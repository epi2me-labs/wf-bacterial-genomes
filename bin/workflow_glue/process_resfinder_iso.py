#!/usr/bin/env python
"""Process resfinder information from output JSON for isolate report sequencing."""


def get_databases(json_file):
    """Get details on database used for resfinder analysis."""
    db_dict = dict()
    databases = json_file["databases"]
    for k, v in databases.items():
        db_dict[k] = {
            "name": v["database_name"],
            "db_version": v["database_version"]
            }
    return db_dict


def get_acquired_data(json_file):
    """Get information on acquired resistace genes found."""
    db_dict = get_databases(json_file)
    acquired_dict = dict()
    acquired = json_file["seq_regions"]
    for k, v in acquired.items():
        gene = k.split(";;")[0]  # For accordion - remove if not using
        if not v["phenotypes"]:
            continue
        acquired_dict[gene] = {
            "gene": v["name"],
            "drugs": v["phenotypes"],
            "database": ", ".join([db_dict[i]["name"] for i in v["ref_database"]]),
            "start": v["query_start_pos"],
            "end": v["query_end_pos"],
            "contig": v["query_id"],
            "identity": round(v["identity"], 3),
            "coverage": round(v["coverage"], 3),
            "pmids": ", ".join(v["pmids"])
            }
    return acquired_dict


def get_point_data(json_file):
    """Get information on point mutations found in resfinder analysis."""
    db_dict = get_databases(json_file)
    point_dict = dict()
    points = json_file["seq_variations"]
    for k, v in points.items():
        gene = k.split(";;")[0]
        if gene not in point_dict:
            point_dict[gene] = list()
        if v["phenotypes"]:
            point_dict[gene].append({
                "gene": k.split(";;")[0],
                "drugs": v["phenotypes"],
                "start": v["ref_start_pos"],
                "end": v["ref_end_pos"],
                "database": db_dict[v["ref_database"]]["name"],
                "aa": v["seq_var"],
                "nuc": f"{v['ref_codon']}>{v['var_codon']}",
                "pmids": ", ".join(v["pmids"])
            })
    return point_dict


def get_drug_class(json_file):
    """Get information on drug classes per drug."""
    pheno_dict = dict()
    phenotypes = json_file["phenotypes"]
    for k, v in phenotypes.items():
        if v["amr_resistant"] and v["amr_species_relevant"]:
            pheno_dict[k] = v["amr_classes"][0]
    return pheno_dict
