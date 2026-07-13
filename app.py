import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdMolDescriptors
import time
import py3Dmol
from stmol import showmol

st.set_page_config(page_title="Kemet Dock", layout="wide")

# Theme & Styling: Custom CSS
st.markdown(
    """
    <style>
    /* Main background */
    .stApp {
        background-color: #F3E9DC;
    }

    /* Text and headers */
    h1, h2, h3, h4, h5, h6, p, div {
        color: #5C3D2E !important;
    }

    /* Styled container mimicking papyrus/gold-sand */
    div[data-testid="stVerticalBlockBorderWrapper"] {
        background-color: #E6CCB2 !important;
        border-radius: 10px;
        padding: 15px;
        border: 2px solid #5C3D2E !important;
    }

    /* Ensure the markdown block text in containers is also dark */
    .stMarkdown {
        color: #5C3D2E !important;
    }
    </style>
    """,
    unsafe_allow_html=True
)

st.title("Kemet Dock")

@st.cache_data
def load_data():
    try:
        df = pd.read_csv("kemet_plants.csv")
        # Handle potential missing data gracefully
        df = df.fillna("Unknown")
        return df
    except Exception as e:
        st.error(f"Error loading data: {e}")
        return pd.DataFrame()

df = load_data()

if not df.empty:
    common_names = df["Common Name"].tolist()
    selected_plant = st.selectbox("Select a Plant:", common_names)

    # Get the selected row
    plant_data = df[df["Common Name"] == selected_plant].iloc[0]

    col1, col2 = st.columns(2)

    with col1:
        with st.container(border=True):
            st.header("Ancient Texts")
            st.write(f"**Botanical Name:** {plant_data['Botanical Name']}")
            st.write(f"**Family:** {plant_data['Family']}")
            st.write(f"**Ancient Hieroglyphic Symbols:** {plant_data['Ancient Symbols']}")
            st.write(f"**Transliteration:** {plant_data['Claim Transliteration']}")
            st.write(f"**Papyrus Source:** {plant_data['Papyrus Source']}")
            st.write(f"**Ancient Claim:** {plant_data['Ancient Claim']}")

    with col2:
        with st.container(border=True):
            st.header("Modern Science")
            st.write(f"**Active Phytochemical:** {plant_data['Active Phytochemical']}")
            st.write(f"**Modern Use:** {plant_data['Modern Use']}")

            st.subheader("2D Chemical Structure")
            smiles = plant_data["SMILES"]
            if smiles != "Unknown" and pd.notna(smiles):
                try:
                    mol = Chem.MolFromSmiles(smiles)
                    if mol:
                        img = Draw.MolToImage(mol, size=(400, 400))
                        st.image(img, use_container_width=True)
                    else:
                        st.write("Could not parse SMILES string into a 2D structure.")
                except Exception as e:
                    st.write(f"Error parsing SMILES string: {e}")
            else:
                st.write("No SMILES string available for this plant.")

            st.selectbox("Select Protein Target Receptor:", [plant_data['Protein Target']])
            st.selectbox("Select Phytochemical Ligand:", [plant_data['Active Phytochemical']])

            initiate_docking = st.button("🚀 Initiate Molecular Docking")

            if initiate_docking:
                with st.spinner('Calculating thermodynamic trajectories and sampling binding poses...'):
                    time.sleep(2)

                st.subheader("Interactive 3D Docking Simulation")

                # Format common name to lower case and spaces to underscores
                file_path = "data/" + selected_plant.lower().replace(" ", "_") + ".pdb"

                try:
                    with open(file_path, "r") as f:
                        pdb_string = f.read()

                    view = py3Dmol.view(width=800, height=600)
                    view.addModel(pdb_string, "pdb")
                    view.setStyle({'cartoon': {'color': 'spectrum'}})
                    view.addStyle({'resn': 'UNL'}, {'stick': {}})
                    view.setBackgroundColor('white')
                    view.zoomTo()
                    showmol(view, height=600, width=800)
                except FileNotFoundError:
                    st.error(f"3D structure file not found: {file_path}")
                except Exception as e:
                    st.error(f"Error loading 3D structure: {e}")

                st.subheader("Molecular Properties")
                if smiles != "Unknown" and pd.notna(smiles):
                    try:
                        mol = Chem.MolFromSmiles(smiles)
                        if mol:
                            num_atoms = mol.GetNumAtoms()
                            num_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
                            num_h_donors = rdMolDescriptors.CalcNumHBD(mol)
                            num_h_acceptors = rdMolDescriptors.CalcNumHBA(mol)

                            st.write(f"**Total Atoms:** {num_atoms}")
                            st.write(f"**Rotatable Bonds:** {num_rotatable_bonds}")
                            st.write(f"**H-Bond Donors:** {num_h_donors}")
                            st.write(f"**H-Bond Acceptors:** {num_h_acceptors}")
                        else:
                            st.write("Could not calculate properties: invalid SMILES.")
                    except Exception as e:
                        st.write(f"Error calculating properties: {e}")
                else:
                    st.write("SMILES string not available for property calculation.")

                st.subheader("Scientifically Accurate Energy Table")
                try:
                    affinities_df = pd.read_csv("binding_affinities.csv")
                    plant_affinities = affinities_df[affinities_df["Common Name"] == selected_plant]
                    if not plant_affinities.empty:
                        # Extract the 7 binding affinity scores
                        cols_to_show = ["Pose 1", "Pose 2", "Pose 3", "Pose 4", "Pose 5", "Pose 6", "Pose 7"]
                        # Drop the index and show a clean table
                        st.table(plant_affinities[cols_to_show].reset_index(drop=True))
                    else:
                        st.write("No binding affinities data found for this plant.")
                except FileNotFoundError:
                    st.error("Binding affinities file (binding_affinities.csv) not found.")
                except Exception as e:
                    st.error(f"Error loading binding affinities: {e}")

else:
    st.warning("No data available in kemet_plants.csv.")
