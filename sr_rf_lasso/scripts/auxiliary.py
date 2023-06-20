
# Auxiliary functions

from functools import reduce
from typing import List

import feyn
import numpy as np
import pandas as pd
from feyn import Model, metrics
from sklearn.preprocessing import MinMaxScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.metrics import accuracy_score, roc_auc_score, precision_recall_curve, f1_score, auc
from sklearn.model_selection import GridSearchCV, train_test_split, StratifiedKFold
from sklearn.pipeline import Pipeline



def agg_sample(df, feature_list, operation):
    """ This function selects columns with each name in feature_list and add a new column based on operation """
    for name in feature_list:
        if operation == "mean":
            df[name] = df.filter(like = name).mean(axis = 1)
            
        elif operation == "std":
            df[name] = df.filter(like = name).std(axis = 1)
            
        elif operation == "median":
            df[name] = df.filter(like = name).median(axis = 1)
            
        elif operation == "sum":
            df[name] = df.filter(like = name).sum(axis = 1)
            
        elif operation == "extreme":
            # Judge which extreme number to take
            Se = np.abs(df.filter(like = name).max(axis = 1)) > np.abs(df.filter(like = name).min(axis = 1))
            array = np.zeros((len(df.index)))
            # Assign max numbers
            array[Se == True] = df.filter(like = name).max(axis = 1)[Se == True]
            # Assign min numbers
            array[Se == False] = df.filter(like = name).min(axis = 1)[Se == False]
            # Add new column to df
            df[name] = array
            
        else:
            print("The input operation term is not included")

            
            
            
            
def modsum(models, train, validation, problem_type):
    """ This function displays the measurement scores of the models """
    
    model_list=[]
    auc_list_train=[]
    auc_list_val=[]
    rmse_train=[]
    rmse_val=[]
    r2_train=[]
    r2_test=[]
    bic_list=[]
    feat_list=[]
    function_list = []
    loss_list=[]
    i=0
    
    if problem_type == 'classification':
        for x in models:
            model_list.append(str(i))
            auc_list_train.append(str(x.roc_auc_score(train).round(2)))
            auc_list_val.append(str(x.roc_auc_score(validation).round(2)))
            bic_list.append(str(x.bic.round(2)))
            feat_list.append(len(x.features))
            function_list.append(str(x.sympify(symbolic_lr=False, symbolic_cat=True, include_weights=False)))
            loss_list.append(x.loss_value)
            i+=1
        df = pd.DataFrame(list(zip(function_list, auc_list_train, auc_list_val, bic_list, feat_list, function_list, loss_list)), 
               columns =['Model', 'AUC Train', 'AUC Test', 'BIC', 'N. Features', 'Functional form', 'Loss'])
    
    elif problem_type == 'regression':
        for x in models:
            model_list.append(str(i))
            rmse_train.append(str(x.rmse(train).round(2)))
            rmse_val.append(str(x.rmse(validation).round(2)))
            r2_train.append(str(x.r2_score(train).round(2)))
            r2_test.append(str(x.r2_score(validation).round(2)))
            bic_list.append(str(x.bic.round(2)))
            feat_list.append(len(x.features))
            function_list.append(str(x.sympify(symbolic_lr=False, symbolic_cat=True, include_weights=False)))
            loss_list.append(x.loss_value)
            i+=1
        df = pd.DataFrame(list(zip(function_list, rmse_train, rmse_val, r2_train, r2_test, 
                                   bic_list, feat_list, function_list, loss_list)), 
               columns =['Model', 'RMSE Train', 'RMSE Test', 'R2 Train', 'R2 Test', 
                         'BIC', 'N. Features', 'Functional form', 'Loss'])
        
    else:
        print('You need to enter a correct problem type, either regression or classification.')

    return(df)




def crossvalidation_as_framework(df, target, n_folds=5, random_state=42, use_sample_weights=True, cv=None,
                                 **kwargs):
    if cv:
        kfold_test = cv
    else:
        kfold_test = StratifiedKFold(n_folds, shuffle=True, random_state=random_state)

    results = ModelResults()

    for i, (train, val) in enumerate(kfold_test.split(df, df[target])):
        train, val = df.iloc[train], df.iloc[val]

        if use_sample_weights:
            sample_weights = np.where(train[target] == 1, np.sum(train[target] == 0)/sum(train[target] == 1), 1)
        else:
            sample_weights = None

        ql = feyn.QLattice(random_seed=42)
        
        # Fit models
        models = ql.auto_run(train, target, sample_weights=sample_weights, **kwargs)


        for j in models:
            results.update(train, val, i, j)

    return results.df





def random_forest_benchmark(data, target, sample_weight=True, param_grid=None, n_folds=5, n_jobs=-1):
    
    # One-hot transformation
    data = pd.get_dummies(data)
    
    # Split data
    train_val, test = train_test_split(data, test_size=0.2, random_state=42, stratify=data[target])
    X_train_val = train_val.iloc[:,:-1]
    X_test = test.iloc[:,:-1]
    y_train_val = train_val[target]
    y_test = test[target]
    
    # Define parameters
    if param_grid:
        param_grid = param_grid
    else:
        param_grid = {'rf__max_depth': [3, 4, 5],
                      'rf__n_estimators': [50, 75, 100],
                      'rf__min_samples_split': [2, 4, 6, 8]}
    
    # Define classifier
    pipeline = Pipeline([('scaler', MinMaxScaler()),
                    ('rf', RandomForestClassifier())])

    grid_search = GridSearchCV(estimator=pipeline, param_grid=param_grid, cv=n_folds, scoring='roc_auc', n_jobs=n_jobs)
    # estimator: classifier
    # n_jobs: number of jobs to run in parallel

    # fit the model with sample weights if needed
    if sample_weight:
        # Assign higher weight to class 1
        sw = np.where(train_val[target] == 1, np.sum(train_val[target] == 0)/sum(train_val[target] == 1), 1)
        grid_search.fit(X_train_val, y_train_val, **{'rf__sample_weight': sw})
    else:
        grid_search.fit(X_train_val, y_train_val)

    
    # Get the best estimator and predict the outcome of test set
    best_estimator = grid_search.best_estimator_
    y_prob = best_estimator.predict_proba(X_test)[:,1]

    return best_estimator, y_prob, y_test




def lasso_benchmark(data, target, sample_weight=True, param_grid=None, n_folds=5, n_jobs=-1):
    # One-hot transformation
    data = pd.get_dummies(data)
    
    # Split data
    train_val, test = train_test_split(data, test_size=0.2, random_state=42, stratify=data[target])
    X_train_val = train_val.iloc[:,:-1]
    X_test = test.iloc[:,:-1]
    y_train_val = train_val[target]
    y_test = test[target]
    
    # Define parameters
    if param_grid:
        param_grid = param_grid
    else:
        param_grid = [1000, 300, 100, 30, 10, 3, 1, .3, .1, .03, .01, .003, .001, .0003, .0001]

    # Define classifier
    pipeline = Pipeline([('scaler', MinMaxScaler()),
                         ('lr', LogisticRegressionCV(Cs=param_grid, penalty='l1', solver='liblinear',
                                                     scoring='roc_auc'))])

    # fit the model with sample weights if needed
    if sample_weight:
        # Assign higher weight to class 1
        sw = np.where(train_val[target] == 1, np.sum(train_val[target] == 0)/sum(train_val[target] == 1), 1)
        pipeline.fit(X_train_val, y_train_val, **{'lr__sample_weight': sw})
    else:
        pipeline.fit(X_train_val, y_train_val)

     # Get the best estimator and predict the outcome of test set
    best_estimator = pipeline.named_steps['lr']
    y_prob = best_estimator.predict_proba(X_test)[:,1]
    return best_estimator, y_prob, y_test




class ModelResults:
    def __init__(self, kind="classification"):
        self.kind = kind

        if self.kind == "classification":
            self.df = pd.DataFrame(columns=['model_structure', 'fold', 'aic', 'bic', 'roc_auc_train',
                                            'accuracy_train', 'roc_auc_val', 'accuracy_val', 'pr_auc', 'f1'])
        elif self.kind == "regression":
            self.df = pd.DataFrame(columns=['model_structure', 'fold', 'aic', 'bic', 'rmse', 'mae', 'r2_score'])
        else:
            raise ValueError("kind must be classification or regression")

    def update(self, train, val, fold, model):
        if model:
            if self.kind == "classification":
                model_precision, model_recall, _ = precision_recall_curve(val[model.output], model.predict(val))
                if model:
                    self.df = self.df.append(pd.DataFrame(data={
                        'model_structure': [str(model.sympify(include_weights=False))],
                        'query_string': [model.to_query_string()],
                        'fold': [fold],
                        'aic': [model.aic],
                        'bic': [model.bic],
                        'roc_auc_train': [model.roc_auc_score(train)],
                        'accuracy_train': [model.accuracy_score(train)],
                        'roc_auc_val': [model.roc_auc_score(val)],
                        'accuracy_val': [model.accuracy_score(val)],
                        'pr_auc': [auc(model_recall, model_precision)],
                        'f1': [f1_score(val[model.output], model.predict(val) > model.accuracy_threshold(train)[0])]
                    }))

            elif self.kind == "regression":
                if model:
                    preds = model.predict(val)
                    self.df = self.df.append(pd.DataFrame(data={
                        'model_structure': [str(model.sympify(include_weights=False))],
                        'query_string': [model.to_query_string()],
                        'fold': [fold],
                        'aic': [model.aic],
                        'bic': [model.bic],
                        'rmse': [feyn.metrics.rmse(val[model.output], preds)],
                        'mae': [feyn.metrics.mae(val[model.output], preds)],
                        'r2_score': [feyn.metrics.r2_score(val[model.output], preds)]
                    }))
                    
