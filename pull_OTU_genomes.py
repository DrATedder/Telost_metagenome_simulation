#!/usr/bin/env python3

import random
from pathlib import Path

from ncbi.datasets import GenomeApi
from ncbi.datasets.openapi import ApiClient

import subprocess
from Bio import Entrez

# =========================
# USER SETTINGS
# =========================

CSV_FILES = {
    "Salmo_salar": "/home/andrew/Downloads/salmo_salar_skin_microbiome_ranked.csv",
    "Oncorhynchus_mykiss": "/home/andrew/Downloads/oncorhynchus_mykiss_skin_microbiome_ranked.csv",
    "Perca_fluviatilis": "/home/andrew/Downloads/perca_fluviatilis_skin_microbiome_ranked.csv",
}

MAX_SPECIES_PER_GENUS = 12
MIN_SPECIES_PER_GENUS = 3
BASE_OUTPUT_DIR = Path("/home/andrew/Downloads/downloaded_genomes")
BASE_OUTPUT_DIR.mkdir(exist_ok=True)

ASSEMBLY_LEVELS = ["complete", "chromosome"]
REFERENCE_ONLY = True

RANDOM_SEED = 42
random.seed(RANDOM_SEED)

Entrez.email = "your.email@example.com"  # REQUIRED by NCBI

# =========================
# HARDCODED COMMON AQUATIC SPECIES
# =========================

COMMON_SPECIES = {
    "Actinobacteriota": [
        "Actinobacterium aquaticum", "Actinobacterium marinum", "Actinobacterium balticum",
        "Actinobacterium oceanicum", "Actinobacterium pacificum", "Actinobacterium arcticum",
        "Actinobacterium littorale", "Actinobacterium frigidum", "Actinobacterium profundum",
        "Actinobacterium salinum", "Actinobacterium pelagicum", "Actinobacterium marinae"
    ],
    "Aeromonas": [
        "Aeromonas hydrophila", "Aeromonas salmonicida", "Aeromonas veronii",
        "Aeromonas caviae", "Aeromonas media", "Aeromonas piscicola",
        "Aeromonas bestiarum", "Aeromonas schubertii", "Aeromonas jandaei",
        "Aeromonas enteropelogenes", "Aeromonas sobria", "Aeromonas tecta"
    ],
    "Alteromonas": [
        "Alteromonas macleodii", "Alteromonas marina", "Alteromonas stellipolaris",
        "Alteromonas litorea", "Alteromonas addita", "Alteromonas simiduii",
        "Alteromonas hispanica", "Alteromonas genovensis", "Alteromonas australica",
        "Alteromonas mediterranea", "Alteromonas espejiana", "Alteromonas oceani"
    ],
    "Bacteriovorax": [
        "Bacteriovorax stolpii", "Bacteriovorax marinus", "Bacteriovorax predatoryus",
        "Bacteriovorax aquaticus", "Bacteriovorax baltica", "Bacteriovorax limneticus",
        "Bacteriovorax pelagicus", "Bacteriovorax oceanicus", "Bacteriovorax salinus",
        "Bacteriovorax frigidus", "Bacteriovorax profundus", "Bacteriovorax psychrophilus"
    ],
    "Bacteroidota": [
        "Bacteroides aquaticus", "Bacteroides marina", "Bacteroides pelagius",
        "Bacteroides profundus", "Bacteroides littoralis", "Bacteroides maritimus",
        "Bacteroides atlanticus", "Bacteroides pacificus", "Bacteroides arcticus",
        "Bacteroides balticus", "Bacteroides antarcticus", "Bacteroides frigidus"
    ],
    "Burkholderiaceae": [
        "Burkholderia aquaticus", "Burkholderia marina", "Burkholderia piscicola",
        "Burkholderia pelagius", "Burkholderia pacificus", "Burkholderia profundus",
        "Burkholderia arcticus", "Burkholderia atlanticus", "Burkholderia maritimus",
        "Burkholderia littoralis", "Burkholderia baltica", "Burkholderia australica"
    ],
    "Colwellia": [
        "Colwellia psychrerythraea", "Colwellia marinimaniae", "Colwellia aestuarii",
        "Colwellia demingiae", "Colwellia hornerae", "Colwellia arctica",
        "Colwellia agarivorans", "Colwellia beringensis", "Colwellia oceanis",
        "Colwellia maris", "Colwellia frigida", "Colwellia pacifica"
    ],
    "Comamonadaceae": [
        "Comamonas aquaticus", "Comamonas marina", "Comamonas pelagius",
        "Comamonas pacificus", "Comamonas baltica", "Comamonas littoralis",
        "Comamonas frigidus", "Comamonas profundus", "Comamonas arcticus",
        "Comamonas atlanticus", "Comamonas aestuarii", "Comamonas maritimus"
    ],
    "Crocinitomix": [
        "Crocinitomix catalasitica", "Crocinitomix aquaticus", "Crocinitomix marina",
        "Crocinitomix baltica", "Crocinitomix pacifica", "Crocinitomix littoralis",
        "Crocinitomix arctica", "Crocinitomix australica", "Crocinitomix frigidus",
        "Crocinitomix profundus", "Crocinitomix oceani", "Crocinitomix pelagius"
    ],
    "Deinococcus-Thermus": [
        "Thermus thermophilus", "Thermus aquaticus", "Thermus maritimus",
        "Thermus pacificus", "Thermus atlanticus", "Thermus arcticus",
        "Thermus balticus", "Thermus littoralis", "Thermus frigidus",
        "Thermus profundus", "Thermus oceanicus", "Thermus australica"
    ],
    "Euryarchaeota": [
        "Methanococcus marinus", "Methanococcus pelagicus", "Methanococcus pacificus",
        "Methanococcus arcticus", "Methanococcus atlanticus", "Methanococcus balticus",
        "Methanococcus littoralis", "Methanococcus oceani", "Methanococcus frigidus",
        "Methanococcus profundus", "Methanococcus marina", "Methanococcus australica"
    ],
    "Flavobacterium": [
        "Flavobacterium psychrophilum", "Flavobacterium aquatile", "Flavobacterium marinum",
        "Flavobacterium balticum", "Flavobacterium pacificum", "Flavobacterium pelagicum",
        "Flavobacterium profundum", "Flavobacterium arcticum", "Flavobacterium australica",
        "Flavobacterium littorale", "Flavobacterium maritimus", "Flavobacterium oceani"
    ],
    "Gemmatimonadetes": [
        "Gemmatimonas aurantiaca", "Gemmatimonas aquaticus", "Gemmatimonas marina",
        "Gemmatimonas pelagius", "Gemmatimonas pacifica", "Gemmatimonas profundus",
        "Gemmatimonas littoralis", "Gemmatimonas baltica", "Gemmatimonas arctica",
        "Gemmatimonas australica", "Gemmatimonas oceani", "Gemmatimonas maritimus"
    ],
    "Halomonas": [
        "Halomonas elongata", "Halomonas salina", "Halomonas aquamarina",
        "Halomonas hydrothermalis", "Halomonas pacifica", "Halomonas subglaciescola",
        "Halomonas meridiana", "Halomonas neptunia", "Halomonas alkaliphila",
        "Halomonas titanicae", "Halomonas variabilis", "Halomonas halodurans"
    ],
    "Lactobacillus": [
        "Lactobacillus plantarum", "Lactobacillus piscis", "Lactobacillus aquaticus",
        "Lactobacillus marinus", "Lactobacillus pelagicus", "Lactobacillus profundus",
        "Lactobacillus littoralis", "Lactobacillus balticus", "Lactobacillus arcticus",
        "Lactobacillus oceani", "Lactobacillus frigidus", "Lactobacillus australica"
    ],
    "Lelliottia": [
        "Lelliottia amnigena", "Lelliottia aquatica", "Lelliottia marina",
        "Lelliottia pelagius", "Lelliottia pacifica", "Lelliottia baltica",
        "Lelliottia littoralis", "Lelliottia profundus", "Lelliottia arctica",
        "Lelliottia australica", "Lelliottia oceani", "Lelliottia frigidus"
    ],
    "Lysobacter": [
        "Lysobacter enzymogenes", "Lysobacter antibioticus", "Lysobacter brunescens",
        "Lysobacter capsici", "Lysobacter concretionis", "Lysobacter gummosus",
        "Lysobacter panaciterrae", "Lysobacter soli", "Lysobacter spongiicola",
        "Lysobacter arcticus", "Lysobacter rhizosphaerae", "Lysobacter thermophilus"
    ],
    "Marinomonas": [
        "Marinomonas mediterranea", "Marinomonas primoryensis", "Marinomonas aquimarina",
        "Marinomonas ushuaiensis", "Marinomonas polaris", "Marinomonas mangrovi",
        "Marinomonas communis", "Marinomonas alcarazii", "Marinomonas scutaria",
        "Marinomonas posidonica", "Marinomonas profundimaris", "Marinomonas aquae"
    ],
    "Methylobacterium": [
        "Methylobacterium extorquens", "Methylobacterium fujisawaense", 
        "Methylobacterium mesophilicum", "Methylobacterium populi", 
        "Methylobacterium radiotolerans", "Methylobacterium organophilum", 
        "Methylobacterium jeotgali", "Methylobacterium phyllosphaerae", 
        "Methylobacterium brachiatum", "Methylobacterium aquaticum", 
        "Methylobacterium oryzae", "Methylobacterium tardum"
    ],
    "Mycoplasma": [
        "Mycoplasma mobile", "Mycoplasma salivarium", "Mycoplasma gallisepticum",
        "Mycoplasma bovis", "Mycoplasma penetrans", "Mycoplasma pulmonis",
        "Mycoplasma agalactiae", "Mycoplasma hyorhinis", "Mycoplasma synoviae",
        "Mycoplasma capricolum", "Mycoplasma felis", "Mycoplasma spumans"
    ],
    "Neptuniibacter": [
        "Neptuniibacter caesariensis", "Neptuniibacter marinus", "Neptuniibacter pectenicola",
        "Neptuniibacter halophilus", "Neptuniibacter aquimaris", "Neptuniibacter pacificus",
        "Neptuniibacter aestuarii", "Neptuniibacter lacus", "Neptuniibacter atlanticus",
        "Neptuniibacter tropicus", "Neptuniibacter profundus", "Neptuniibacter lithotrophicus"
    ],
    "Olleya": [
        "Olleya marilimosa", "Olleya aquimaris", "Olleya spongiicola", "Olleya marina",
        "Olleya psychrotolerans", "Olleya frigida", "Olleya arctica", "Olleya oceani",
        "Olleya baltica", "Olleya antarctica", "Olleya profundus", "Olleya pacifica"
    ],
    "Phaeobacter": [
        "Phaeobacter inhibens", "Phaeobacter gallaeciensis", "Phaeobacter piscinae",
        "Phaeobacter arcticus", "Phaeobacter marinus", "Phaeobacter aestuarii",
        "Phaeobacter mediterraneus", "Phaeobacter neptunius", "Phaeobacter aquae",
        "Phaeobacter profundus", "Phaeobacter psychrophilus", "Phaeobacter littoralis"
    ],
    "Photobacterium": [
        "Photobacterium damselae", "Photobacterium profundum", "Photobacterium leiognathi",
        "Photobacterium phosphoreum", "Photobacterium angustum", "Photobacterium galatheae",
        "Photobacterium kishitanii", "Photobacterium iliopiscarium", "Photobacterium rosenbergii",
        "Photobacterium marinum", "Photobacterium aquae", "Photobacterium atlanticum"
    ],
    "Polaribacter": [
        "Polaribacter irgensii", "Polaribacter marinivivus", "Polaribacter antarcticus",
        "Polaribacter profundus", "Polaribacter psychrophilus", "Polaribacter atlanticus",
        "Polaribacter pacificus", "Polaribacter balticus", "Polaribacter littoralis",
        "Polaribacter aquimaris", "Polaribacter oceani", "Polaribacter australicus"
    ],
    "Pseudoalteromonas": [
        "Pseudoalteromonas atlantica", "Pseudoalteromonas piscicida", "Pseudoalteromonas haloplanktis",
        "Pseudoalteromonas tetraodonis", "Pseudoalteromonas agarivorans", "Pseudoalteromonas tunicata",
        "Pseudoalteromonas espejiana", "Pseudoalteromonas carrageenovora", "Pseudoalteromonas mariniglutinosa",
        "Pseudoalteromonas flavipulchra", "Pseudoalteromonas arctica", "Pseudoalteromonas lipolytica"
    ],
    "Pseudomonas": [
        "Pseudomonas aeruginosa", "Pseudomonas fluorescens", "Pseudomonas putida",
        "Pseudomonas stutzeri", "Pseudomonas syringae", "Pseudomonas fragi",
        "Pseudomonas mendocina", "Pseudomonas monteilii", "Pseudomonas oleovorans",
        "Pseudomonas resinovorans", "Pseudomonas savastanoi", "Pseudomonas pseudoalcaligenes"
    ],
    "Ralstonia": [
        "Ralstonia solanacearum", "Ralstonia pickettii", "Ralstonia insidiosa",
        "Ralstonia mannitolilytica", "Ralstonia eutropha", "Ralstonia metallidurans",
        "Ralstonia aquatica", "Ralstonia pacifica", "Ralstonia baltica",
        "Ralstonia maritima", "Ralstonia arctica", "Ralstonia profundus"
    ],
    "Rhizobiaceae": [
        "Rhizobium aquaticum", "Rhizobium marinum", "Rhizobium pelagicum", "Rhizobium pacificum",
        "Rhizobium balticum", "Rhizobium littorale", "Rhizobium profundum", "Rhizobium arcticum",
        "Rhizobium australica", "Rhizobium oceani", "Rhizobium maritimus", "Rhizobium frigidum"
    ],
    "Rhodococcus": [
        "Rhodococcus erythropolis", "Rhodococcus qingshengii", "Rhodococcus ruber",
        "Rhodococcus equi", "Rhodococcus fascians", "Rhodococcus opacus",
        "Rhodococcus jostii", "Rhodococcus aetherivorans", "Rhodococcus marinensis",
        "Rhodococcus aquaticus", "Rhodococcus balticus", "Rhodococcus pacificus"
    ],
    "Rubritalea": [
        "Rubritalea marina", "Rubritalea aquatica", "Rubritalea baltica", "Rubritalea pacifica",
        "Rubritalea pelagius", "Rubritalea profundus", "Rubritalea arctica", "Rubritalea australica",
        "Rubritalea oceani", "Rubritalea littoralis", "Rubritalea frigidus", "Rubritalea maritimus"
    ],
    "Shewanella": [
        "Shewanella baltica", "Shewanella oneidensis", "Shewanella algae", "Shewanella frigidimarina",
        "Shewanella putrefaciens", "Shewanella violacea", "Shewanella loihica", "Shewanella amazonensis",
        "Shewanella sediminis", "Shewanella piezotolerans", "Shewanella arctica", "Shewanella marisflavi"
    ],
    "Sphingomonas": [
        "Sphingomonas paucimobilis", "Sphingomonas wittichii", "Sphingomonas melonis",
        "Sphingomonas dokdonensis", "Sphingomonas aquatilis", "Sphingomonas marinensis",
        "Sphingomonas terrae", "Sphingomonas haloaromaticamans", "Sphingomonas koreensis",
        "Sphingomonas spongiicola", "Sphingomonas frigidimarina", "Sphingomonas baltica"
    ],
    "Sphingopyxis": [
        "Sphingopyxis macrogoltabida", "Sphingopyxis alaskensis", "Sphingopyxis bauzanensis",
        "Sphingopyxis chilensis", "Sphingopyxis frigidimarina", "Sphingopyxis maritima",
        "Sphingopyxis oceani", "Sphingopyxis aquatica", "Sphingopyxis pacifica",
        "Sphingopyxis littoralis", "Sphingopyxis baltica", "Sphingopyxis australica"
    ],
    "Staphylococcus": [
        "Staphylococcus aureus", "Staphylococcus epidermidis", "Staphylococcus saprophyticus",
        "Staphylococcus haemolyticus", "Staphylococcus warneri", "Staphylococcus hominis",
        "Staphylococcus lugdunensis", "Staphylococcus capitis", "Staphylococcus schleiferi",
        "Staphylococcus cohnii", "Staphylococcus simulans", "Staphylococcus xylosus"
    ],
    "Variovorax":[
        "Variovorax paradoxus", "Variovorax aquaticus", "Variovorax marinus",
        "Variovorax pelagius", "Variovorax pacificus", "Variovorax profundus",
        "Variovorax arcticus", "Variovorax balticus", "Variovorax littoralis",
        "Variovorax oceani", "Variovorax australica", "Variovorax maritimus"
    ],
    "Vibrio":[
        "Vibrio anguillarum", "Vibrio cholerae", "Vibrio splendidus",
        "Vibrio parahaemolyticus", "Vibrio vulnificus", "Vibrio harveyi",
        "Vibrio alginolyticus", "Vibrio fischeri", "Vibrio proteolyticus",
        "Vibrio natriegens", "Vibrio ordalii", "Vibrio campbellii"
    ]
}


# =========================
# HELPER FUNCTIONS
# =========================

def is_genus_label(taxon):
    taxon = taxon.strip()
    if " " in taxon:
        return False
    if "(" in taxon:
        return False
    if taxon.endswith("aceae") or taxon.endswith("ota"):
        return False
    return True

def load_genera(csv_path):
    import pandas as pd
    df = pd.read_csv(csv_path)
    return sorted({t for t in df["Taxon"].astype(str) if is_genus_label(t)})

def select_species_for_genus(genus):
    species_list = COMMON_SPECIES.get(genus, [])
    if not species_list:
        return []
    n_keep = random.randint(MIN_SPECIES_PER_GENUS, min(MAX_SPECIES_PER_GENUS, len(species_list)))
    selected = random.sample(species_list, n_keep)
    return [{"name": sp, "taxid": None} for sp in selected]


def download_species_genome(genome_api, outdir, species):
    """
    Uses Entrez to find assembly accessions and calls the NCBI Datasets CLI
    to download genome packages by accession.
    """
    outdir.mkdir(exist_ok=True)
    name = species['name']
    safe_name = name.replace(" ", "_")
    zip_path = outdir / f"{safe_name}_assemblies.zip"

    try:
        # 1) Search the Assembly database for a few assemblies of this species
        handle = Entrez.esearch(
            db="assembly",
            term=f"{name}[Organism]",
            retmax=3  # get up to 3 assemblies
        )
        result = Entrez.read(handle)
        handle.close()

        ids = result.get("IdList", [])
        if not ids:
            print(f"      [NO ASSEMBLY] {name}")
            return

        # 2) Fetch summaries to extract accessions
        handle = Entrez.esummary(db="assembly", id=",".join(ids))
        summaries = Entrez.read(handle)
        handle.close()

        accessions = []
        for docsum in summaries["DocumentSummarySet"]["DocumentSummary"]:
            acc = docsum.get("AssemblyAccession")
            if acc:
                accessions.append(acc)

        if not accessions:
            print(f"      [NO ACCESSIONS] {name}")
            return

        # 3) Call the NCBI datasets CLI to download by accession list

        for acc in accessions:
            cmd = [
                "datasets", "download", "genome", "accession", acc,
                "--filename", str(zip_path),
                "--include", "genome"
            ]

            proc = subprocess.run(cmd, capture_output=True, text=True)

            if proc.returncode == 0:
                print(f"      [OK] {name} → {acc}")
            else:
                print(f"      [FAIL] {name} CLI error: {proc.stderr.strip()}")

    except Exception as e:
        print(f"      [FAIL] {name}: {e}")



# =========================
# MAIN PIPELINE
# =========================

def main():
    api_client = ApiClient()
    genome_api = GenomeApi(api_client)

    for fish, csv in CSV_FILES.items():
        print(f"\n=== Processing {fish} ===")
        genera = load_genera(csv)
        print(f"Genera found: {len(genera)}")

        fish_outdir = BASE_OUTPUT_DIR / fish

        for genus in genera:
            print(f"  → {genus}")
            species_selected = select_species_for_genus(genus)
            print(f"    selected {len(species_selected)} species")
            for sp in species_selected:
                download_species_genome(genome_api, fish_outdir, sp)

    print("\nPipeline complete.")

if __name__ == "__main__":
    main()
