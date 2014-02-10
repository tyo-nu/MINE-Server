#!/usr/bin/python
#
# Call the CDMI API to get data needed to run the autorecon algorithm
# Needed data includes:
#
# For likelihood computations:
# - List of substems and the sequences of their members
# - List of OTU organisms
# - Gene neighborhoods for OTU organisms
#
# For optimization:
# - Reactions [only want mass/charge-balanced ones]
# - Metabolites
# -[Growth data? ]

# The CDMI_API is for "well-trodden paths" functions
# CDMI_EntityAPI is for ER functions (all_entities_..., get_Relationship_....)
from biokbase.cdmi.client import CDMI_API, CDMI_EntityAPI
from urllib2 import URLError, HTTPError
import urllib
import sys
import operator #for itemgetter

try:
    import json
except ImportError:
    sys.path.append('simplejson-2.3.3')
    import simplejson as json

def getFieldFromEntity(seedEntity, fieldName):
    ''' Get field "fieldName" from the entity seedEntity, the result of a 
    get_entity_xxx or an all_entities_xxx call to the CDM. The function returns
    a list of entries in field fieldName.

    The entities are dictionaries from ID to a dictionary of key-value pairs
    where keys are whatever you tell it you want. '''

    if seedEntity is None:
        sys.stderr.write("INTERNAL ERROR: Provided seedEntity was None - usually this means you were searching for something that doesnt exist in the database\n")
        raise ValueError
    # Check for an error I seem to make all the time and yell at me in a USEFUL way
    if not isinstance(seedEntity, dict):
        sys.stderr.write("INTERNAL ERROR: getFieldFromEntity expects a dictionary - perhaps you meant to call getFieldFromRelationship?\n")
        raise ValueError
    f = []
    for entry in seedEntity:
        if fieldName not in seedEntity[entry]:
            sys.stderr.write("INTERNAL ERROR: Field name %s not found in provided entity\n" %(fieldName))
            raise ValueError
        f.append(seedEntity[entry][fieldName])
    return f


def getFieldFromRelationship(seedRelationship, fieldName, objtype):
    '''
    INPUTS:
    seedRelationship: The result of one of the get_relationship_xxx functions
    fieldName: The field you want to extract from the object.
    objtype: "TO", "REL", or "FROM"

    OUTPUTS:
    A list (in the same order as the list from the get_relationship function)
    of the values with the specified field name.

    The get_relationship_xxx functions return lists of lists.
    The first list is a list of all the links
    The second list has three dictionaries in it: the TO dictionary, the REL dictionary and the FROM dictionary
    describing properties on either end of the link and of the link itself...

    If you want to  maintain duplicicate relationships (many-to-one, one-to-many, many-to-many), this function should be called at
    least twice (once on each end of the relationship, or once on an end and once in the middle)..
    '''

    if seedRelationship is None:
        sys.stderr.write("INTERNAL ERROR: The provided relationship was None - usually this means you were searching for something that doesn't exist in the database.\n")
        raise ValueError
    objidx = None
    if objtype.lower() == "from":
        objidx = 0
    elif objtype.lower() == "rel":
        objidx = 1
    elif objtype.lower() == "to":
        objidx = 2
    else:
        sys.stderr.write("INTERNAL ERROR: In getFieldFromRelationship - objtype must be TO, REL, or FROM\n")
        raise ValueError
    if not isinstance(seedRelationship, list):
        sys.stderr.write("INTERNAL ERROR: getFieldFromRelationship expects a list - perhaps you meant to call getFieldFromEntity?\n")
        raise ValueError
    # Unravel
    f = []
    for entry in seedRelationship:
        # TO CHECK: Is it possible to have one of the links lead to nothing?
        # Check field name validity - if it links to something there has to be the data request there
        # or else something is wrong.
        if fieldName not in entry[objidx]:
            sys.stderr.write("INTERNAL ERROR: Field name %s not found in provided relationship\n" %(fieldName))
            raise ValueError
        f.append(entry[objidx][fieldName])
    return f

def subsystemFids(count, config):
    '''Query the CDMI for a list of FIDs in the subsystems'''

    cdmi = CDMI_API(config["cdmi_url"])
    cdmi_entity = CDMI_EntityAPI(config["cdmi_url"])

    # Get the genes that are in subsystems and in OTUs.
    ssdict = dict()
    start = 0
    done = False
    while not done:
        subdict = cdmi_entity.all_entities_Subsystem(start, count, ["id"])
        ssdict.update(subdict)
        start += count
        if len(subdict) < count:
            done = True
    ssids = getFieldFromEntity(ssdict, "id")

    # Now lets get a list of FIDs within those subsystems
    # Break the complete list into smaller sub-lists to avoid timeouts
    start = 0
    increment = 8
    end = start + increment - 1
    counter = len(ssids)
    ssfids = []
    while counter > 0:
        try:
            ssfiddict = cdmi.subsystems_to_fids(ssids[start:end], [])
        except HTTPError as e:
            if increment > 1:
                increment = increment / 2
                end = start + increment - 1
            sys.stderr.write("caught '%s' error, increment is now %d\n" %(e.reason, increment))
            continue
        for key in ssfiddict:
            for ssfid in ssfiddict[key]:
                ls = ssfiddict[key][ssfid]
                for arr in ls:
                    if len(arr) > 1:
                        gl = arr[1]
                        for l in gl:
                            ssfids.append(l)
                            
        # Move to next sub-list
        start += increment
        end += increment
        if end >= len(ssids):
            end = len(ssids) - 1
        counter -= increment

    # Uniquify!
    return list(set(ssfids))

def getDlitFids(count, config):
    '''Get a list of FIDs with direct literature evidence (dlits) from the CDM'''

    cdmi = CDMI_API(config["cdmi_url"])
    cdmi_entity = CDMI_EntityAPI(config["cdmi_url"])
    pubdict = dict()
    start = 0
    done = False
    while not done:
        subdict = cdmi_entity.all_entities_Publication(start, count, ["id"])
        pubdict.update(subdict)
        start += count
        if len(subdict) < count:
            done = True

    pubids = getFieldFromEntity(pubdict, "id")
    sys.stderr.write("%d pubids\n" %(len(pubids)))
    pub2seq = cdmi_entity.get_relationship_Concerns(pubids, [], [], ["id"])
    pubseqs = getFieldFromRelationship(pub2seq, "id", "to")
    sys.stderr.write("%d pubseqs\n" %(len(pubseqs)))
    seq2fids = cdmi_entity.get_relationship_IsProteinFor(pubseqs, [], [], ["id"])
    fids = getFieldFromRelationship(seq2fids, "id", "to")
    return fids

def filterFidsByOtus(fidlist, otus, config):
    '''
    Obsolete (I think this isn't used any more)

    Given a list of representative organism IDs (OTUs) and a list of
    FIDs, returns only those FIDs found in an OTU.'''

    cdmi_entity = CDMI_EntityAPI(config["cdmi_url"])

    # Identify the organism belonging to each fid
    # If this fails to find an organism we don't want it anyway...
    orgdict = cdmi_entity.get_relationship_IsOwnedBy(fidlist, [], [], ["id"])
    flist = getFieldFromRelationship(orgdict, "from_link", "rel")
    olist = getFieldFromRelationship(orgdict, "id", "to")

    fids = []
    for ii in range(len(olist)):
        if olist[ii] in otus:
            fids.append(flist[ii])
    return fids

def filterFidsByOtusBetter(fidsToRoles, rolesToFids, oturepsToMembers, config):
    '''Attempt to do a more intelligent filtering of FIDs by OTU.

    Given all FIDs attached to a role in the unfiltered set we do the following:
    
    Initialize KEEP
    For each OTU and each role:
       If role is found in the representative, add to KEEP and continue;
       Otherwise, iterate over other genomes.
           If role is found in one other genome, add to KEEP and continue;

    This process should make our calculation less sensitive to the choice of OTUs...

    '''

    cdmi_entity = CDMI_EntityAPI(config["cdmi_url"])

    # Identify the organism belonging to each fid
    # If this fails to find an organism we don't want it anyway...
    fidlist = fidsToRoles.keys()
    orgdict = []
     # Break the complete list into smaller sub-lists to avoid timeouts
    start = 0
    increment = 5000
    end = start + increment - 1
    counter = len(fidlist)
    while counter > 0:
        try:
            od = cdmi_entity.get_relationship_IsOwnedBy(fidlist[start:end], [], [], ["id"])
        except HTTPError as e:
            if increment > 1:
                increment = increment / 2
                end = start + increment - 1
            sys.stderr.write("caught '%s' error, increment is now %d\n" %(e.reason, increment))
            continue
        orgdict.extend(od)
        start += increment
        end += increment
        if end >= len(fidlist):
            end = len(fidlist) - 1
        counter -= increment
    fidlist = getFieldFromRelationship(orgdict, "from_link", "rel")
    orglist = getFieldFromRelationship(orgdict, "id", "to")
    fidToOrg = {}
    for ii in range(len(fidlist)):
        fidToOrg[fidlist[ii]] = orglist[ii]
    
    keptFidsToRoles = {}
    keptRolesToFids = {}
    # If the OTUs are comprehensive this should be empty.
    missingRoles = []

    # For each OTU
    for oturep in oturepsToMembers:
        # for each role
        for role in rolesToFids:
            fidlist = rolesToFids[role]
            keepFid = None
            keepRole = None
            for fid in fidlist:
                # This can happen due to MOL issues
                if fid not in fidToOrg:
                    continue
                org = fidToOrg[fid]
                # If the organism is the representative we keep it and go to the next role
                if org == oturep:
                    keepFid = fid
                    keepRole = role
                    break
                # Otherwise look at the rest of the list (note that I just pick one without really paying
                # attention to WHICH one...). We save them in case there are no examples of the role in the
                # representative organism, but continue on anyway.
                if org in oturepsToMembers[oturep]:
                    keepFid = fid
                    keepRole = role
            if keepFid is not None:
                if keepFid in keptFidsToRoles:
                    keptFidsToRoles[keepFid].append(keepRole)
                else:
                    keptFidsToRoles[keepFid] = [ keepRole ]
                if keepRole in keptRolesToFids:
                    keptRolesToFids[keepRole].append(keepFid)
                else:
                    keptRolesToFids[keepRole] = [ keepFid ]

    missingRoles = list(set(rolesToFids.keys()) - set(keptRolesToFids.keys()))

#    print oturepsToMembers
#    print missingRoles
#    print keptRolesToFids

    return keptFidsToRoles, keptRolesToFids, missingRoles

def getOtuGenomeDictionary(count, config):
    '''Obtain a dictionary from OTU representatives to all genomes in the OTU'''
    cdmi = CDMI_API(config["cdmi_url"])
    # Get list of OTUs
    otulist = getOtuGenomeIds(count, config)
    otudict = cdmi.otu_members(otulist[0])
    return otudict

def fidsToRoles(fidlist, config):
    '''Given a list of FIDs return a dictionary from FID to the list of roles the encoding gene
    performs and a dictionary from roles to the FIDs performing them'''

    cdmi = CDMI_API(config["cdmi_url"])
    cdmi_entity = CDMI_EntityAPI(config["cdmi_url"])
    
    # Break the complete list into smaller sub-lists to avoid timeouts
    start = 0
    increment = 1000
    end = start + increment - 1
    counter = len(fidlist)
    fidsToRoles = {}
    rolesToFids = {}
    while counter > 0:
        try:
            roledict = cdmi_entity.get_relationship_HasFunctional(fidlist[start:end], [], [], ["id"])
        except HTTPError as e:
            if increment > 1:
                increment = increment / 2
                end = start + increment - 1
            sys.stderr.write("caught '%s' error, increment is now %d\n" %(e.reason, increment))
            continue
        flist = getFieldFromRelationship(roledict, "from_link", "rel")
        rolelist = getFieldFromRelationship(roledict, "id", "to")
        for ii in range(len(flist)):
            # We have to use sets here because a bug(?) in get_relationship_HasFunctional allows multiple identical
            # links between fids and roles.
            # See for example what happens when you call it on g.9647.peg.2332
            if flist[ii] in fidsToRoles:
                fidsToRoles[flist[ii]].add(rolelist[ii])
            else:
                fidsToRoles[flist[ii]] = set([rolelist[ii]])
            if rolelist[ii] in rolesToFids:
                rolesToFids[rolelist[ii]].add(flist[ii])
            else:
                rolesToFids[rolelist[ii]] = set([flist[ii]])
                
        # Move to next sub-list
        start += increment
        end += increment
        if end >= len(fidlist):
            end = len(fidlist) - 1
        counter -= increment
        
    # Convert back to lists to not break other functions.
    for f in fidsToRoles:
        fidsToRoles[f] = list(fidsToRoles[f])
    for r in rolesToFids:
        rolesToFids[r] = list(rolesToFids[r])
    return fidsToRoles, rolesToFids

def fidsToSequences(fidlist, config):
    '''Given a list of FIDs, returns a dictionary from FID to its amino acid sequence.
    Features with no amino acid sequence are discarded.'''
    cdmi = CDMI_API(config["cdmi_url"])
    fidlist = list(set(fidlist))
    start = 0
    increment = 5000
    end = start + increment - 1
    counter = len(fidlist)
    seqs = {}
    while counter > 0:
        try:
            ps = cdmi.fids_to_protein_sequences(fidlist[start:end])
        except HTTPError as e:
            if increment > 1:
                increment = increment / 2
                end = start + increment - 1
            sys.stderr.write("caught '%s' error, increment is now %d\n" %(e.reason, increment))
            continue
        seqs.update(ps)
        
        # Move to next sub-list
        start += increment
        end += increment
        if end >= len(fidlist):
            end = len(fidlist) - 1
        counter -= increment
    
    return seqs

def genomesToPegs(genomes, config):
    '''Given a list fogenome IDs, returns a list of FIDs for protein-encoding genes in the specified genomes'''

    cdmi_entity = CDMI_EntityAPI(config["cdmi_url"])
    fiddict = cdmi_entity.get_relationship_IsOwnerOf(genomes, [], [], ["id", "feature_type"])
    fidlist = getFieldFromRelationship(fiddict, "id", "to")
    typelist = getFieldFromRelationship(fiddict, "feature_type", "to")
    # We want protein-encoding genes only (not e.g. operons, operators, etc...)
    # The type of protein-encoding genes is CDS now but will possibly be changed to peg later...
    pegs = []
    for ii in range(len(fidlist)):
        if typelist[ii] == "peg" or typelist[ii] == "CDS":
            pegs.append(fidlist[ii])    
    return pegs

def getOtuGenomeIds(count, config):
    '''Query the CDMI for a list of OTU genomes (returns a list of OTUs and a list of only
    prokaryote OTUs)'''

    # Get the complete list of OTUs.
    cdmi_entity = CDMI_EntityAPI(config["cdmi_url"])
    otudict = dict()
    start = 0
    done = False
    while not done:
        subdict = cdmi_entity.all_entities_OTU(start, count, ["id"])
        otudict.update(subdict)
        start += count
        if len(subdict) < count:
            done = True

    # Find out if a OTU is marked as representative and if it is prokaryotic.
    otuids = getFieldFromEntity(otudict, "id")
    gendict = cdmi_entity.get_relationship_IsCollectionOf(otuids, [], ["representative"], ["id", "prokaryotic"])
    isrep = getFieldFromRelationship(gendict, "representative", "rel")
    isprok = getFieldFromRelationship(gendict, "prokaryotic", "to")
    genomeid = getFieldFromRelationship(gendict, "id", "to")
    prokotus = []
    otus = []
    for ii in range(len(genomeid)):
        if int(isrep[ii]) == 1 and int(isprok[ii]) == 1:
            prokotus.append(genomeid[ii])
        if int(isrep[ii]) == 1:
            otus.append(genomeid[ii])
    return otus, prokotus

################
# Input: genomes - List of genome IDs
# Output 1: Tuples of gene neighborhood information
#   (contig_id, feature_id, start_location, strand)
#   for each input OTU
#   The output is sorted by contig id, then by start location
#
# Output 2: A dictionary from fid to roles
################
def getGenomeNeighborhoodsAndRoles(genomes, config):
    cdmi_entity = CDMI_EntityAPI(config["cdmi_url"])

    pegs = genomesToPegs(genomes)
    # Get contigs
    fidlocdict = cdmi_entity.get_relationship_IsLocatedIn(pegs, [], ["begin", "dir"], ["id"])
    fids = getFieldFromRelationship(fidlocdict, "from_link", "rel")
    begins = getFieldFromRelationship(fidlocdict, "begin", "rel")
    dirs = getFieldFromRelationship(fidlocdict, "dir", "rel")
    cids = getFieldFromRelationship(fidlocdict, "id", "to")

    tuplist = []
    for ii in range(len(cids)):
        tuplist.append( (cids[ii], fids[ii], int(begins[ii]), dirs[ii]) )
    # Sort by contig first, then by start location.
    tuplist = sorted(tuplist, key=operator.itemgetter(0,2))

    # Now lets get the role for all of these IDs
    # Note that a single protein can have multiple roles.
    roledict = cdmi_entity.get_relationship_HasFunctional(fids, [], [], ["id"])
    fids = getFieldFromRelationship(roledict, "from_link", "rel")
    roles = getFieldFromRelationship(roledict, "id", "to")
    fidToRoles = {}
    rolesToFids = {}
    for ii in range(len(fids)):
        if fids[ii] in fidToRoles:
            fidToRoles[fids[ii]].append(roles[ii])
        else:
            fidToRoles[fids[ii]] = [ roles[ii] ]
        if roles[ii] in rolesToFids:
            rolesToFids[roles[ii]].append(fids[ii])
        else:
            rolesToFids[roles[ii]] = [ fids[ii] ]
    return tuplist, fidToRoles

def complexRoleLinks(count, config):
    '''
    OBSOLETE - will be replaced by Chris's roles_to_reactions() function

    Query the CDM for a list of links from complexes to roles. Returns a dictionary from
    role to a list of complexes and a dictionary from complexes to a list of roles.

    Only roles listed as "required" are included in the links.'''

    # Get a list of complexes
    cdmi_entity = CDMI_EntityAPI(config["cdmi_url"])
    cplxdict = dict()
    start = 0
    done = False
    while not done:
        subdict = cdmi_entity.all_entities_Complex(start, count, ["id"])
        cplxdict.update(subdict)
        start += count
        if len(subdict) < count:
            done = True
    cplxlist = getFieldFromEntity(cplxdict, "id")

    # Get a list of roles linked to those complexes
    roledict = cdmi_entity.get_relationship_IsTriggeredBy(cplxlist, [], ["optional"], ["id"])
    cplx = getFieldFromRelationship(roledict, "from_link", "rel")
    opt = getFieldFromRelationship(roledict, "optional", "rel")
    role = getFieldFromRelationship(roledict, "id", "to")
    complexToRequiredRoles = {}
    requiredRolesToComplex = {}
    for ii in range(len(cplx)):
        # For now - we don't want to deal with the "optional" components. I'm not sure how I'd incorporate them into a likelihood calculation anyway.
        if int(opt[ii]) == 1:
            continue
        # Note - this becomes an all-AND GPR - (role1 AND role2 AND ... )
        if cplx[ii] in complexToRequiredRoles:
            complexToRequiredRoles[cplx[ii]].append(role[ii])
        else:
            complexToRequiredRoles[cplx[ii]] = [ role[ii] ]
        if role[ii] in requiredRolesToComplex:
            requiredRolesToComplex[role[ii]].append(cplx[ii])
        else:
            requiredRolesToComplex[role[ii]] = [ cplx[ii] ]
    return complexToRequiredRoles, requiredRolesToComplex

def reactionComplexLinks(count, config):
    '''
    OBSOLETE - will be replaced by Chris's roles_to_reactions() function
    
    Query the CDM for a list of links from reactions to complexes. Returns a dictionary
    from reactions to lists of complexes performing them and from complexes to lists of reactions
    they perform.
    '''

    cdmi_entity = CDMI_EntityAPI(config["cdmi_url"])

    # The API was recently changed to use model IDs and to not use the reactions_to_complexes
    # but use the ER model instead.
    # I reflect that here...
    rxndict = dict()
    start = 0
    done = False
    while not done:
        subdict = cdmi_entity.all_entities_Reaction(count, ["id"])
        rxndict.update(subdict)
        start += count
        if len(subdict) < count:
            done = True
    rxns = getFieldFromEntity(rxndict, "id")
    cplxdict = cdmi_entity.get_relationship_IsStepOf(rxns, [], [], ["id"])
    rxnlist = getFieldFromRelationship(cplxdict, "from_link", "rel")
    cplxlist = getFieldFromRelationship(cplxdict, "id", "to")
    
    rxnToComplex = {}
    complexToRxn = {}
    for ii in range(len(rxnlist)):
        if rxnlist[ii] in rxnToComplex:
            rxnToComplex[rxnlist[ii]].append(cplxlist[ii])
        else:
            rxnToComplex[rxnlist[ii]] = [ cplxlist[ii] ]
        if cplxlist[ii] in complexToRxn:
            complexToRxn[cplxlist[ii]].append(rxnlist[ii])
        else:
            complexToRxn[cplxlist[ii]] = [ rxnlist[ii] ]

    return rxnToComplex, complexToRxn
