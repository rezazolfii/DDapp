import pandas as pd
import numpy as np
from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
from rdkit.Chem import AllChem
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit import DataStructs
from rdkit.Chem import rdFingerprintGenerator
from padelpy import padeldescriptor
from sklearn.preprocessing import StandardScaler
from sklearn.manifold import TSNE
from sklearn.cluster import MiniBatchKMeans
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import os
import subprocess




class GeneratingFingerprints:
    def __init__(self):
        pass

    def fetch_user_target(self, target_name):
        target = new_client.target
        target_query = target.search(target_name)
        self.targets = pd.DataFrame.from_dict(target_query)

        if self.targets.empty:
            raise ValueError("No targets found with the specified name.")
        self.targets = self.targets.loc[self.targets["organism"] == "Homo sapiens"]

        return self.targets


    def extract_user_target(self, user_selected_index):

        df_target_proteins = None

        for n in user_selected_index:

            # Retrieve the selected target's ChEMBL ID
            selected_target = self.targets.loc[n, 'target_chembl_id']

            # Fetch activity data
            activity = new_client.activity
            res = activity.filter(target_chembl_id=selected_target).filter(standard_type="IC50")

            # Convert the activity results to a DataFrame
            df_target_protein = pd.DataFrame.from_dict(res)

            if df_target_protein.empty:
                raise ValueError(f"No activity data found for index {n}.")
                continue

            if df_target_proteins is None:
                df_target_proteins = df_target_protein
            else:
                df_target_proteins = pd.concat([df_target_proteins, df_target_protein])

        self.df_target_proteins = df_target_proteins.reset_index().drop(columns = "index")

        return self.df_target_proteins


    def extracting_fingerprints(self):

        selection = ['molecule_chembl_id', 'canonical_smiles']
        df_selection = self.df_target_proteins[selection].dropna().reset_index().drop(columns = "index")

        def lipinski(smiles, verbose=False):
            """imputs a dataframe of smile strings to return a dataframe with lipiski descriptors"""

            moldata = []
            for elem in smiles:
                mol = Chem.MolFromSmiles(elem)
                if mol is not None:
                    moldata.append(mol)
                elif verbose:
                    print(f"Invalid SMILES: {elem}")
            baseData = np.empty((0, 4))

            for mol in moldata:
                desc_MolWt = Descriptors.MolWt(mol)
                desc_MolLogP = Descriptors.MolLogP(mol)
                desc_NumHDonors = Lipinski.NumHDonors(mol)
                desc_NumHAcceptors = Lipinski.NumHAcceptors(mol)

                row = np.array([desc_MolWt, desc_MolLogP, desc_NumHDonors, desc_NumHAcceptors])

                baseData = np.vstack([baseData, row])

            columnNames = ["MW", "LogP", "NumHDonors", "NumHAcceptors"]
            descriptors = pd.DataFrame(data=baseData, columns=columnNames)
            return descriptors

        df_lipinski = lipinski(df_selection["canonical_smiles"])

        def morgen(smiles):

            mols = [Chem.MolFromSmiles(i) for i in smiles]
            calc = MoleculeDescriptors.MolecularDescriptorCalculator([x[0] for x in Descriptors._descList])
            desc_names = calc.GetDescriptorNames()
            Mol_descriptors = []

            for mol in mols:
                mol = Chem.AddHs(mol)
                descriptors = calc.CalcDescriptors(mol)
                Mol_descriptors.append(descriptors)

            df_morgen = pd.DataFrame(Mol_descriptors, columns=desc_names)

            return df_morgen

        df_morgen = morgen(df_selection['canonical_smiles'])

        def desc_calc():
            # Performs the descriptor calculation
            bashCommand = "java -Xms2G -Xmx2G -Djava.awt.headless=true -jar ./PaDEL-Descriptor/PaDEL-Descriptor.jar -removesalt -standardizenitro -fingerprints -descriptortypes ./PaDEL-Descriptor/PubchemFingerprinter.xml -dir ./ -file descriptors_output.csv"
            process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
            output, error = process.communicate()
            os.remove('molecule.smi')

        df_selection['canonical_smiles'].to_csv('smiles_input.smi', index=False, header=False)

        smiles_input_path = 'smiles_input.smi'
        output_csv = 'descriptors_output.csv'

        padeldescriptor(mol_dir=smiles_input_path, d_file=output_csv, fingerprints=True)
        df_padel = pd.read_csv(output_csv)

        self.df_final = pd.concat([df_selection["molecule_chembl_id"], df_lipinski, df_morgen, df_padel], axis = 1).drop(columns = "Name")

        return self.df_final


    def fit_and_predict(self):

        ids = self.df_final["molecule_chembl_id"]

        X = self.df_final.drop(columns = "molecule_chembl_id")

        invalid = ['BCUT2D_MWHI', 'BCUT2D_MWLOW', 'BCUT2D_CHGHI', 'BCUT2D_CHGLO', 'BCUT2D_LOGPHI', 'BCUT2D_LOGPLOW', 'BCUT2D_MRHI', 'BCUT2D_MRLOW']
        if [i for i, n in X.isna().sum().to_dict().items() if n != 0] == invalid:
            X = X.drop(columns = invalid)
        else:
            raise KeyError("Too bad, try again")

        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)

        n_samples = X_scaled.shape[0]
        perplexity = min(30, n_samples - 1)

        tsne = TSNE(n_components=2, perplexity=perplexity)
        X_reduced = tsne.fit_transform(X_scaled)

        def elbow_method(X, min_cluster, max_clusters, step_size):

            wcss = []
            clusters = list(range(min_cluster, max_clusters + 1, step_size))

            for i in clusters:
                kmeans = MiniBatchKMeans(n_clusters = i, random_state=100)
                kmeans.fit(X_scaled)
                wcss.append(kmeans.inertia_)

            inertia_diff = np.diff(np.diff(wcss)) # Compute the differences between each inertia value

            optimal_clusters = np.argmin(inertia_diff) + 2  # Adding 2 because we lose two points with np.diff twice

            return optimal_clusters

        if X_reduced.shape[0] <= 40:
            min_cluster = 1
            step_size = 1
        elif X_reduced.shape[0] > 40 and X_reduced.shape[0] <= 200:
            min_cluster = 5
            step_size = 2
        elif X_reduced.shape[0] > 200 and X_reduced.shape[0] <= 1200:
            min_cluster = 20
            step_size = 5
        else:
            min_cluster = 30
            step_size = 10

        optimal_clusters = elbow_method(X_reduced, min_cluster, X_reduced.shape[0], step_size)

        clustering = KMeans(n_clusters = optimal_clusters, random_state=100)
        predicted_clusters = clustering.fit_predict(X_reduced)

        y = pd.concat([ids, pd.DataFrame(predicted_clusters, columns=["predicted_cluster"])], axis = 1)

        fig, ax = plt.subplots()
        scatter = ax.scatter(X_reduced[:, 0], X_reduced[:, 1], c=y["predicted_cluster"])
        ax.set_title("t-SNE with KMeans")

        return y, fig, optimal_clusters  # Return the plot and number of clusters
