import torch
import torch.nn as nn
import torch.optim as optim
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from torch.utils.data import DataLoader, Dataset
from sklearn.metrics import roc_auc_score, f1_score, accuracy_score, average_precision_score, precision_score




class MLP(nn.Module):
    """Define the MLP network"""
    def __init__(self, input_size, hidden_size, output_size, l2_reg):
        super(MLP, self).__init__()
        self.fc1 = nn.Linear(input_size, hidden_size)
        self.fc2 = nn.Linear(hidden_size, hidden_size) 
        self.fc3 = nn.Linear(hidden_size, output_size)
        self.init_weights()
        self.relu = nn.ReLU()
        self.l2_reg = l2_reg
        
    def forward(self, x):
        x = self.fc1(x)
        x = self.relu(x)
        x = self.fc2(x)
        x = self.relu(x)
        x = self.fc3(x)
        return x
    
    def init_weights(self):
        for param in self.parameters():
            if param.dim() > 1:
                nn.init.xavier_uniform_(param)
            else:
                nn.init.constant_(param, 0)
    
    def l2_loss(self):
        l2_loss = 0.0
        for param in self.parameters():
            l2_loss += torch.sum(torch.square(param))
        return l2_loss
    
    
    
# Define a custom dataset class for loading tabular data
class CustomDataset(Dataset):
    """Transform numpy array to pytorch tensor"""
    def __init__(self, data, labels):
        # Transform the data from numpy version to pytorch tensor
        self.data = torch.from_numpy(data.values).float()
        self.labels = torch.from_numpy(labels.values).long()
        
    def __len__(self):
        return len(self.data)
    
    def __getitem__(self, idx):
        return self.data[idx], self.labels[idx]
    
    
    
    
# Define the training function
def train(model, dataloader, criterion, optimizer, device, bi_class=True):
    """This function trains the MLP model"""
    model.train() # Set up training mode
    running_loss = 0.0 # Initialize loss
    y_train_true = []
    y_train_scores = []
    
    for inputs, labels in dataloader: # Iterations
        inputs, labels = inputs.to(device), labels.to(device)
        optimizer.zero_grad()
        outputs = model(inputs) # =model.forward(inputs)
        softmax = nn.Softmax(dim=1)
        if bi_class:
            scores = softmax(outputs)[:, 1]
        else:
            scores = softmax(outputs)
        y_train_true.extend(labels.cpu().numpy())
        y_train_scores.extend(scores.detach().numpy())
        loss = criterion(outputs, labels)
        l2_loss = model.l2_loss() * model.l2_reg
        loss_sum = loss + l2_loss
        loss_sum.backward() # Backpropogation
        optimizer.step() # Update parameters
        
        # Update the running loss and predictions
        running_loss += loss.item()
        
    # Calculate overall loss and AUC-ROC of an epoch
    epoch_loss = running_loss / (len(dataloader.dataset)/inputs.size(0)) 
    if bi_class:
        auc_score = roc_auc_score(y_train_true, y_train_scores)
    else:
        auc_score = roc_auc_score(y_train_true, y_train_scores, average='weighted', multi_class='ovr')
    
    return epoch_loss, auc_score




# Define the evaluation function for binary classification
def evaluate(model, dataloader, criterion, device, bi_class=True):
    """This function is for model prediction and evaluation"""
    model.eval() # Set model to evaluation mode
    y_true = []
    y_probs = []
    y_preds = []
    running_loss = 0
    
    with torch.no_grad():
        for inputs, labels in dataloader:
            inputs, labels = inputs.to(device), labels.to(device)
            outputs = model(inputs)
            softmax = nn.Softmax(dim=1)
            if bi_class:
                scores = softmax(outputs)[:, 1]
            else:
                scores = softmax(outputs)
                # Assign the class with the largest probability
                pred = torch.argmax(scores, dim=1)
                y_preds.extend(pred.detach().numpy())
                
            loss = criterion(outputs, labels)
            
            # Update the running loss and predictions
            running_loss += loss.item()
            
            y_true.extend(labels.cpu().numpy())
            y_probs.extend(scores.detach().numpy())
            
    
    # Validation loss
    epoch_loss = running_loss / (len(dataloader.dataset)/inputs.size(0)) 
    
    if bi_class:     
        # Calculate AUC-ROC
        auc_roc = roc_auc_score(y_true, y_probs)

        # Calculate AP (average precision)
        AP = average_precision_score(y_true, y_probs)

        # Set different thresholds
        thresholds = np.linspace(0, 1, 101)
        f1_scores = []
        for threshold in thresholds:
            y_pred = (y_probs >= threshold).astype(int)
            f1_scores.append(f1)

        # Set the threshold to the value that produces the largest accuracy for prediction
        y_pred = (y_probs >= thresholds[f1_scores.index(max(f1_scores))]).astype(bool)

        # Calculate ACC and F1
        accuracy = accuracy_score(y_true, y_pred)
        f1 = f1_score(y_true, y_pred)
        T = thresholds[f1_scores.index(max(f1_scores))]
        
    else:
        auc_roc = roc_auc_score(y_true, y_probs, average='weighted', multi_class='ovr')
        
        # Assign the class with the largest probability
        scores = torch.argmax(scores, dim=1)
        accuracy = accuracy_score(y_true, y_preds)
        f1 = f1_score(y_true, y_preds, average='weighted', zero_division=0)
        precision = precision_score(y_true, y_preds, average='weighted', zero_division=0)
        
    if bi_class:
        return auc_roc, AP, accuracy, f1, T, epoch_loss
    else:
        return auc_roc, accuracy, f1, precision, epoch_loss





def conditions(df):
    """This function transforms categorical target to integers"""
    if df['logFC'] == 'Non_sign':
        return 0
    elif float(df['logFC']) < 0:
        return 1
    else:
        return 2

    

    
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
            
            

       
    

def compute_saliency(model, x):
    """This function computes the gradients of the output with respect to the input features"""
    # Set x as the variable (the parameters are constant)
    x = x.clone().detach().requires_grad_(True)

    # Prediction using the saved model
    class_probs = model(x)
    
    # Compute the derivative of the resulting class with respect to x (How the result is changing with the change of x)
    # 1) Obtain class 
    target = torch.tensor([class_probs.argmax()])
    # 2) Extract class probability score
    target_class_prob = torch.gather(class_probs, 0, target)
    # 3) Get class gradient
    target_class_prob.backward()
    
    # Calculate feature importance of x
    saliency = x.grad.abs().numpy()
    
    # Normalize between 0 and 1, and the sum is 1
    saliency = saliency / sum(saliency)
    
    return saliency





def visualize_saliency(saliency, feature_names, gene_name, data):
    """This function visualizes the saliency map for a given input instance"""
    indices = np.argsort(saliency)[::-1]
    plt.bar(range(data.iloc[:,:-1].shape[1]), saliency[indices], color=['green'], align="center")
    plt.title("{}{}".format("MLP Saliency Map of ", gene_name), fontsize=12)
    plt.xlabel("Feature Name", fontsize=12)
    plt.xticks(range(data.iloc[:,:-1].shape[1]), data.iloc[:,:-1].columns[indices], rotation=90)
    plt.ylabel("Saliency Score", fontsize=12)
    plt.show()
    
    
    
def saliency_heatmap(model, data, class_index, bi_class=True, onlyDEG=False):
    """This function uses the trained model to predict genes, and produces the saliency heatmap of truly-predicted genes"""
    # Prepare the true data points
    X = torch.tensor(data.iloc[:, :-1].values, dtype=torch.float).requires_grad_(True)
    y = torch.tensor(data.iloc[:, -1].values, dtype=torch.long)
    y_pred = model(X)
    softmax = nn.Softmax(dim=1)
    y_pred = softmax(y_pred)
    y_pred = torch.argmax(y_pred, dim=1)
    index = np.where((y_pred == class_index) & (y_pred == y))
    X = X[index]
    y = y[index]

    # Collect the corresponding gene names
    gene_names = []
    for i in range(len(index[0])):
        gene_name = data.index[index[0][i]]
        gene_names.append(gene_name)
    
    if len(gene_names) != 0:
        # Create heatmap for all the genes
        score_matrix = np.zeros((len(gene_names), len(data.iloc[:,:-1].columns)))
        for i in range(len(gene_names)):
            saliency = compute_saliency(model, X[i])
            gene = gene_names[i]
            score_matrix[i,:] = saliency           
        feature_names = data.columns[:-1]
        
        heatmap = plt.imshow(score_matrix, cmap='coolwarm', aspect='auto')
        colorbar = plt.colorbar(heatmap)
        colorbar.set_label('Normalized saliency score')
        plt.xticks(range(data.iloc[:,:-1].shape[1]), feature_names, rotation=90)
        if bi_class:
            if onlyDEG:
                if class_index == 0:
                    plt.ylabel("True DR-DEGs")
                else:
                    plt.ylabel("True UR-DEGs")
            else:
                if class_index == 0:
                    plt.ylabel("True NonDEGs")
                else:
                    plt.ylabel("True DEGs")
        else:
            if class_index == 0:
                plt.ylabel("True NonDEGs")
            elif class_index == 1:
                plt.ylabel("True DR-DEGs")
            elif class_index == 2:
                plt.ylabel("True UR-DEGs")
            else:
                raise ValueError("class_index can only be 0 or 1 or 2.")
        plt.savefig('../figures/saliency_map.png', dpi=300, bbox_inches="tight")
#         plt.show()
        return score_matrix, gene_names
    else:
        print("No truely predicted genes in this class.")