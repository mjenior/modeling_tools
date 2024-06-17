#!/usr/bin/env python

import numpy
import pickle
from Bio import GenBank, SeqIO
from copy import deepcopy
from cobra.flux_analysis import pfba
import cobra, cobra.io, cobra.manipulation



def gbk2faa(gbk_file):

    all_entries = []
    with open(gbk_file, 'r') as GBFile:

        GBcds = SeqIO.InsdcIO.GenBankCdsFeatureIterator(GBFile)

        for cds in GBcds:
            if cds.seq is not None:
                cds.id = cds.name

                if cds.annotations['product'] is not None:
                    cds.description = '_'.join(cds.annotations['product'].split()) + '|' + cds.annotations['genename']
                else:
                    cds.description = 'unknown|' + cds.annotations['genename']

                all_entries.append(cds)

    # write file
    SeqIO.write(all_entries, '{}.faa'.format(gbk_file[:-3]), 'fasta')

    return gbk_file.rstrip('gbk') + 'faa'


# Find exchange reaction ID for a metabolite
def _find_exchange(model, metID):

    for rxn in model.metabolites.get_by_id(metID).reactions:
        if len(rxn.metabolites) == 1:
            return rxn.id
        

# Scale each active exchange back and examine its influence on objective flux
def _find_primary_sources(model, fraction=0.1):

    c_sources = []
    #n_sources = []
    solution = model.optimize()
    objVal = solution.objective_value
    
    # Parse exchange flux samples for imported metabolites
    for rxn in model.boundary:
        flux = solution.fluxes[rxn.id]
        if flux >= -1e-6: continue # Skip exported byproducts or unused reactions

        # Test for disproportionate effect on objective
        old_bound = rxn.lower_bound
        rxn.lower_bound = flux * fraction 
        new_objVal = model.slim_optimize(error_value=0)
        rxn.lower_bound = old_bound # Reset bound
        
        # Calculate the degree of change to objective value
        corrected_flux = abs(flux) / (new_objVal / objVal) if new_objVal != objVal else 1
    
        # Normalize elemental component contributions
        metabolite = rxn.reactants[0]
        for element in metabolite.elements.keys():
            element_supply = float(metabolite.elements[element]) * corrected_flux
            if element_supply > 0.: element_supply = numpy.log(element_supply)
                
            # Identify largest sources of main elements
            if element == 'C' and element_supply > 0.0:
                c_sources.append([metabolite.id, element_supply])
            #elif element == 'N' and element_supply > 0.0:
            #    n_sources.append([metabolite.id, element_supply])
    
    # Rank by largest contributions
    def getKey(item): return item[1]
    c_sources = sorted(c_sources, reverse=True, key=getKey)
    primary_C = c_sources[0][0]
    #n_sources = sorted(n_sources, reverse=True, key=getKey)
    #primary_N = n_sources[0][0]
    
    return primary_C


def _calc_min_flux(model, min_yield):

    model.objective = 'Biomass'
    product_fluxes = pfba(model).fluxes

    substrateID = _find_primary_sources(model)
    substrateExch = _find_exchange(model, substrateID)
    mw_substrate = model.metabolites.get_by_id(substrateID).formula_weight

    yield_mol = -product_fluxes['Biomass'] / product_fluxes[substrateExch]
    yield_g = (yield_mol * 1000.0) / mw_substrate # Mass yields
    yield_ratio = (min_yield / yield_g) if yield_g > min_yield else min_yield

    objFlux = model.slim_optimize()
    min_objFlux = yield_ratio * objFlux
    
    return min_objFlux


def create_strain_model(model, geneDict, strainID, threshold=0.5):

    print('Building strain-specific model for', strainID, 'using the reference model', model.name)
    strain_model = deepcopy(model)
    strain_model.id = strainID
    strain_model.name = strainID
    strain_flux_min = _calc_min_flux(strain_model, threshold)
    print(round(strain_flux_min, 4), 'optimal objective flux minimum\n')

    # Read and parse BLAST results, save discordant hits
    no_hit_model = []
    renameDict = {}
    for gene in strain_model.genes:
        try:
            renameDict[gene.id] = geneDict[gene.id]
        except:
            no_hit_model.append(gene)
    cobra.manipulation.modify.rename_genes(strain_model, renameDict)
    strain_model.repair()

    no_hit_ratio = round((len(no_hit_model) / len(strain_model.genes)) * 100., 2)
    no_hit_ratio = ' (' + str(no_hit_ratio) + '%)\n'
    print(len(no_hit_model), 'model genes did not have an ortholog in reference species model', no_hit_ratio)

    # Identify and save those genes needed for growth
    print('Checking fluxes of non-ortholog reactions...')
    removed_genes = []
    removed_rxns = []
    saved_rxns = []
    disassociated_genes = []
    for gene in no_hit_model:
            
        for rxn in strain_model.genes.get_by_id(gene.id).reactions:
            rxn.bounds = (0.,0.)
            
        test = strain_model.slim_optimize(error_value=0.0)
        if test >= strain_flux_min:
            removed_genes.append(gene.id)
            for rxn in gene.reactions: removed_rxns.append(rxn.id)
            cobra.manipulation.delete.remove_genes(strain_model, [gene], remove_reactions=True)
        else:
            for rxn in gene.reactions:
                saved_rxns.append(rxn.id)
                rxn.bounds = model.reactions.get_by_id(rxn.id).bounds # reset bounds
                rxn.annotation['strain_conversion'] = 'lacks gene ortholog, required for near optimal growth'
            disassociated_genes.append(gene.id)
            cobra.manipulation.delete.remove_genes(strain_model, [gene], remove_reactions=False)

    strain_model.repair()
    print(len(disassociated_genes), 'non-orthologous genes needed for growth minimum, removing GPRs for', len(saved_rxns), 'reactions')
    strain_model.annotation['GPRs removed'] = ', '.join(saved_rxns)
    strain_model.annotation['Gene annotations removed'] = ', '.join(disassociated_genes)

    print(len(removed_genes), 'remaining non-orthologous genes with', len(removed_rxns), 'reactions removed')
    strain_model.annotation['Genes + reactions removed'] = ', '.join(removed_genes)

    # Check for growth
    print('\nOriginal model objective flux:', round(model.slim_optimize(),5))
    print('Strain-specific model objective flux:', round(strain_model.slim_optimize(),5))
    change = round(model.slim_optimize() - strain_model.slim_optimize(),5)
    if change < 0:
        print('\tINCREASED objective flux by:', abs(change), '\n')
    elif change > 0:
        print('\tDECREASED objective flux by:', change, '\n')
    else:
        print('\tNo change in objective flux\n')

    print('New model has', len(strain_model.genes), 'total genes, scaffold model has', len(model.genes))
    print('New model has', len(strain_model.reactions), 'total reactions, scaffold model has', len(model.reactions))
    print('New model has', len(strain_model.metabolites), 'total metabolites, scaffold model has', len(model.metabolites),'\n')

    return strain_model


def gene_conversion(blast_out_1, blast_out_2, in_fasta, out_fasta='check_orfs.faa'):
    
    found_hit = set()
    conversionDict1 = {}
    with open(blast_out_1,'r') as blast:
        for line in blast:
            line = line.split()
            queryID = line[0].split('|')[0]
            found_hit |= set([queryID])
            refID = line[1].split('|')[0]
            conversionDict1[refID] = queryID
    
    conversionDict2 = {}
    with open(blast_out_2,'r') as blast:
        for line in blast:
            line = line.split()
            queryID = line[0].split('|')[0]
            refID = line[1].split('|')[0]
            found_hit |= set([refID])
            conversionDict2[queryID] = refID

    checkORFs = []
    putative = 0
    with open(in_fasta,'r') as orfs:
        with open(out_fasta,'w') as outFile:
            for line in orfs:
                if line[0] == '>':
                    info = line.split('|')
                    orf = info[0].replace('>','')
                    annotation = info[-1].strip()
                    if not orf in found_hit:
                        checkORFs.append(orf)
                        seq = orfs.readline()
                        outFile.write(line)
                        outFile.write(seq + '\n')
                        if annotation in ['hypothetical_protein','putative_protein']: 
                            putative += 1

    for x in conversionDict2.keys():
        try:
            test = conversionDict1[x]
        except:
            conversionDict1[x] = conversionDict2[x]

    print(len(set(checkORFs)), 'total ORFs had no hits across references (including', putative, 'putative/hypothetical ORFs)\n\tFailed matches written to:', out_fasta, '\n')
    conversionDict1['check_orfs'] = checkORFs

    return conversionDict1


def fill_from_universal(model, gbkFile, conversionDict):

    print('Loading reference files...')
    #with open("/Users/mjenior/Desktop/seed_to_bigg/ec_to_bigg.reactions.pkl", "rb") as p:
    with open("/home/databases/ec_to_bigg.reactions.pkl", "rb") as p:
        EC_BiGG_Dict = pickle.load(p)
    #universal = cobra.io.load_json_model('/Users/mjenior/Desktop/seed_to_bigg/universal_model_cobrapy.json')
    universal = cobra.io.load_json_model('/home/databases/universal_model_cobrapy.json')

    print('Parsing strain-specific genes...\n')
    checkORFs = set(conversionDict['check_orfs'])

    parser = GenBank.RecordParser()
    record = parser.parse(open(gbkFile))
    
    locus = 'start'
    failed = []
    updated = []
    added = []
    for feature in record.features:
        for qualifier in feature.qualifiers:
            if qualifier.key == "/locus_tag=":
                locus = qualifier.value.replace('\"','')
            elif qualifier.key == "/product=":
                annotation = qualifier.value.replace('\"','')
                
            if locus in checkORFs and qualifier.key == "/EC_number=": 
                ec_number = qualifier.value.replace('\"','')
                try:
                    rxn_id = EC_BiGG_Dict[ec_number]
                except:
                    continue
                
                try:
                    if model.reactions.get_by_id(rxn_id).gene_reaction_rule == "":
                        model.reactions.get_by_id(rxn_id).gene_reaction_rule = locus
                    else:
                        model.reactions.get_by_id(rxn_id).gene_reaction_rule += f" or {locus}"   
                    updated.append(f'{locus} --> {rxn_id}')
                except:
                    try:
                        model.add_reactions([deepcopy(universal.reactions.get_by_id(rxn_id))])
                        model.reactions.get_by_id(rxn_id).gene_reaction_rule = locus
                        model.genes.get_by_id(locus).name = annotation
                        added.append(f'{locus} --> {rxn_id}')
                    except:
                        failed.append(f'{locus} --> {rxn_id}')
    
    if len(added) >= 1:
        print('Created new reactions with GPRs:')
        for x in set(added): print('\t', x)
        print('\n')
    if len(updated) >= 1:
        print('Added to GPRs for already existing reactions:')
        for y in set(updated): print('\t', y)
        print('\n')
    if len(failed) >= 1:
        print('Failed to find reactions for these BiGG IDs:')
        for z in set(failed): print('\t', z)
        print('\n')
    
    print('Model now has', len(model.genes), 'total genes')
    print('Model now has', len(model.reactions), 'total reactions')
    print('Model now has', len(model.metabolites), 'total metabolites\n')

    return model


# Checks remaining ORFs against KEGG
def integrate_unique_genes(model, fasta):

    print('\nChecking remaining failed annotations against KEGG (can take a while)...')
    zoolander_path = "/home/databases"
    #zoolander_path = "/Users/mjenior/Desktop/conversion_refs" # local
    with open(f"{zoolander_path}/gene_to_bigg.pkl", "rb") as p: kegg_bigg_dict = pickle.load(p)
    universal = cobra.io.load_json_model(f"{zoolander_path}/universal_model_cobrapy.json")
    
    if not "tsv" in fasta:
        kegg_hits = _run_blast(fasta, f"{zoolander_path}/KEGG_prokaryotes_pep")
    else:
        kegg_hits = fasta

    hits = {}
    with open(kegg_hits, "r") as inFile:
        for line in inFile:
            line = line.split()
            ko = line[1].split("|")[1]
            if ko != "none":
                locus = line[0].split("|")[0]
                gene = line[1].split("|")[0]
                annotation = line[1].split("|")[2]
                hits[locus] = [ko, gene, annotation]

    added = []
    for locus in hits.keys():
        try:
            annotation = hits[locus]
            rxnID = kegg_bigg_dict[annotation[0]]
        except:
            continue
        
        try:
            test = model.reactions.get_by_id(rxnID)
        except:
            try:
                model.add_reactions([deepcopy(universal.reactions.get_by_id(rxnID))])
                model.reactions.get_by_id(rxnID).gene_reaction_rule = locus
                model.genes.get_by_id(locus).name = annotation[1]
                added.append(f'{locus} --> {rxnID} | {annotation[2]}')
            except:
                continue

    if len(added) >= 1:
        print('\nIdentified new reactions and genes:')
        for x in set(added): print('\t', x)
    else:
        print('\nNo additional metabolic reactions found.')

    return model
