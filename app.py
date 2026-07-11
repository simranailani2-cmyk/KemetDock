import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw

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
            st.write(f"**Protein Target:** {plant_data['Protein Target']}")
            st.write(f"**PDB ID:** {plant_data['PDB ID']}")

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
else:
    st.warning("No data available in kemet_plants.csv.")
