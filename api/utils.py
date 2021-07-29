from flask import current_app as app
from rdkit.Chem import AllChem
from rdkit.Chem.AllChem import MolFromSmiles, MolToInchiKey


def get_smiles_from_mol_string(mol_string):
    """Convert a molfile in string format to a SMILES string."""
    mol = AllChem.MolFromMolBlock(mol_string)
    smiles = AllChem.MolToSmiles(mol)

    return smiles


def get_extra_info(db, core_db, ref_db, compounds):
    """Look up compounds in core database to get spectra, fingerprints,
    and DB Links."""

    cpd_ids = []
    final_compounds = {}
    for cpd in compounds:
        cpd_ids.append(cpd['_id'])

        mol = MolFromSmiles(cpd['SMILES'])
        inchi_key = MolToInchiKey(mol)
        key_prefix = inchi_key.split('-')[0]

        ref_cpds = ref_db.data.find({'Inchikey': {'$regex': r'^' + key_prefix}})
        best_ref_cpd = _get_best_ref_cpd(ref_cpds)
        if best_ref_cpd:
            cpd['Cross_References'], cpd['Names'] = _get_xrefs(best_ref_cpd)
        else:
            cpd['Cross_References'] = {}
            cpd['Names'] = []

        final_compounds[cpd['_id']] = cpd

    core_compounds = core_db.compounds.find({'_id': {'$in': cpd_ids}})

    for core_cpd in core_compounds:
        cpd_id = core_cpd['_id']
        final_cpd = final_compounds[cpd_id]
        final_cpd['Mass'] = core_cpd['Mass']
        final_cpd['Charge'] = core_cpd['Charge']
        final_cpd['Formula'] = core_cpd['Formula']
        final_cpd['Inchikey'] = core_cpd['Inchikey']
        final_cpd['logP'] = core_cpd['logP']
        final_cpd['RDKit_fp'] = core_cpd['RDKit_fp']
        final_cpd['len_RDKit_fp'] = core_cpd['len_RDKit_fp']
        final_cpd['Spectra'] = core_cpd['Spectra']
        final_cpd['MINE_id'] = core_cpd['MINE_id']
        final_cpd['KEGG_id'] = core_cpd['KEGG_id']

    final_compounds = [final_compounds[cpd_id] for cpd_id in final_compounds]

    return final_compounds


def _get_best_ref_cpd(ref_cpds):
    """Select reference compound with most cross references."""
    best_ref_cpd = None
    pubchem_id = None
    for ref_cpd in ref_cpds:
        if 'pubchem_id' in ref_cpd:
            pubchem_id = ref_cpd['pubchem_id']

        if 'cross_references' in ref_cpd:
            if not best_ref_cpd:
                best_ref_cpd = ref_cpd
            else:
                if len(ref_cpd['cross_references']) > len(best_ref_cpd['cross_references']):
                    best_ref_cpd = ref_cpd

    if best_ref_cpd and 'pubchem_id' not in best_ref_cpd and pubchem_id:
        best_ref_cpd['pubchem_id'] = pubchem_id

    return best_ref_cpd


def _get_xrefs(cpd_dict):
    """Get cross-references and names for compound."""
    if not cpd_dict:
        return None, None
    xrefs = {}
    names = set()
    if 'cross_references' in cpd_dict:
        for xref in cpd_dict['cross_references']:
            if _description_is_valid(xref['description']):
                xrefs[xref['source']] = xref['source_id']
                for name in xref['description'].split('||'):
                    if name.lower() not in [n.lower() for n in names]:
                        names.add(name)

    if 'pubchem_id' in cpd_dict:
        xrefs['pubchem_id'] = str(cpd_dict['pubchem_id'])

    return (xrefs, list(names))


def _description_is_valid(description):
    """Check if cross-reference is valid by its description."""
    invalid_descriptions = [
        'secondary/obsolete/fantasy identifier',
        None,
        '',
    ]
    return description not in invalid_descriptions


def extract_enzyme_pathway_from_kegg_data(data):
    """Get enzyme and pathway info from KEGG api response for compound."""
    enzymes = []
    pathways = []
    enzymes_active = False
    pathways_active = False
    for line in data.split('\n'):
        if line.startswith('PATHWAY'):
            pathways_active = True
        elif line.startswith('ENZYME'):
            pathways_active = False
            enzymes_active = True
        elif line.startswith('DBLINKS') or line.startswith('BRITE'):
            enzymes_active = False

        if enzymes_active:
            for enz in line.split():
                enz = enz.strip()
                if enz != 'ENZYME':
                    enzymes.append(enz)

        if pathways_active:
            pathway = line.split('map')[-1].strip()
            pathways.append(pathway)

    return enzymes, pathways
