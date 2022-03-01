## Import packages required for the search
import dill
import math
import numpy
import re
import pandas

## Function to parse json or xml files from the API responses
def getParsedOutput(apiLink, respType):
    import requests
    import json
    resp = requests.get(apiLink)
    searchResult = resp.text
    if respType == 'json':
        # Parse the json text to a dictionary object
        parsedJson = json.loads(searchResult)
        return(parsedJson)
    if respType == 'xml':
        import xmltodict
        parsedXml = xmltodict.parse(searchResult)
        return(parsedXml)
    else:
        print("Only works with JSON and XML formats for now")


######################## Get API results ########################
# In this section we
# - search the APIs of the chosen repositories with particular search terms, with or without filters,
# - count the number of results and
# - save the information to be reviewed for relevance to the initial search
###### Array express ######
# Free text search
arrayexpressLinkAll = 'https://www.ebi.ac.uk/arrayexpress/json/v3/experiments?keywords=huntington+blood+human'
arrayexpressJsonAll = getParsedOutput(apiLink=arrayexpressLinkAll, respType='json')
arrayexpressCountAll = arrayexpressJsonAll['experiments']['total'] # 8 results
# Search with filters, free text for terms without filters
arrayexpressLinkFilt = 'https://www.ebi.ac.uk/arrayexpress/json/v3/experiments?keywords=huntington+blood&species=%22homo%20sapiens%22'
arrayexpressJsonFilt = getParsedOutput(apiLink=arrayexpressLinkFilt, respType='json')
arrayexpressCountFilt = arrayexpressJsonFilt['experiments']['total'] # 8 results
# Experiment details to review
# Relevant details are in 'accession', 'name', 'organism', 'description' and 'samplecharacteristic'
for i in range(arrayexpressCountAll):
    aeExptId = arrayexpressJsonAll['experiments']['experiment'][i]['accession']
    aeExptName = arrayexpressJsonAll['experiments']['experiment'][i]['name']
    aeExptOrganism = arrayexpressJsonAll['experiments']['experiment'][i]['organism'][0]
    aeExptDescription = arrayexpressJsonAll['experiments']['experiment'][i]['description'][0]['text']
    if i == 0:
        aeDF = pandas.DataFrame({'Repo':'ArrayExpress','Id':aeExptId, 'Name':aeExptName, 'Organism':aeExptOrganism, 'Description':aeExptDescription}, index=[0])
    else:
        aeDF = aeDF.append({'Repo':'ArrayExpress','Id':aeExptId, 'Name':aeExptName, 'Organism':aeExptOrganism, 'Description':aeExptDescription}, ignore_index=True)
aeDF.to_csv('arrayExpressAllTerms.csv', index=False)
# Take in information from 'samplecharacteristic' manually because each experiment has different values in this location
categoryValuesDisease = arrayexpressJsonAll['experiments']['experiment'][1]['samplecharacteristic'][4]['value'] # 4 disease, 9 organism part

###### DbGaP ######
# Uses the Entrez API to serve information like most NCBI resources
# Free text search
gapSearchLinkAll = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gap&term=Huntington+blood+human'
gapSearchResultAll = getParsedOutput(gapSearchLinkAll, 'xml')
gapSearchCountAll = gapSearchResultAll['eSearchResult']['Count'] # 21 results
# Search with filters, free text for terms without filters
gapSearchLinkFilt = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gap&term=Huntington[disease]+blood+human'
gapSearchResultFilt = getParsedOutput(gapSearchLinkFilt, 'xml')
gapSearchCountFilt = gapSearchResultFilt['eSearchResult']['Count'] # 1 results
# The Entrez API retrieves the ids of datasets that match the search and a separate query must be made to get
# experimental details of those ids
# Get all ids that matched the search
gapSearchAllIdsLink = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gap&retmax=' + gapSearchCountAll + '&term=Huntington+blood+human'
gapSearchAllIdsResult = getParsedOutput(gapSearchAllIdsLink, 'xml')
# Make a loop here with a retmax of 50 ids per query
gapSraIds = gapSearchAllIdsResult['eSearchResult']['IdList']['Id']
gapNumSplits = math.ceil(len(gapSraIds) / 50)
gapSraSplitIds = numpy.array_split(gapSraIds, gapNumSplits)
gapSraSummary = []
for i in range(gapNumSplits):
    tmpIds = gapSraSplitIds[i]
    tmpGapSummaryLink = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=sra&id=' + ','.join(tmpIds)
    tmpGapSummaryResult = getParsedOutput(tmpGapSummaryLink, 'xml')
    for j in range(len(gapSraSplitIds[i])):
        gapSraSummary.append(tmpGapSummaryResult['eSummaryResult']['DocSum'][j]['Item'][0]['#text'])
# Experiment details to review
# The XML output from the Entrez API is not parsable so here we use regex to retrieve the relevant details for review
for i in range(len(gapSraSummary)):
    status = re.sub('^.* status="(.*)".*$', '\\1', gapSraSummary[i])
    if status != 'withdrawn':
        gapExptId = re.sub('^.*\\<Study acc="(.*)" name="[A-Z][a-z].*$', '\\1', gapSraSummary[i])
        gapExptName = re.sub('^.*\\<Study acc="(.*)" name="(.*)"\\/\\>\\<Organism taxid.*$', '\\2', gapSraSummary[i])
        gapExptOrganism = re.sub('^.*ScientificName="(.*)"\\/\\>\\<Sample.*$', '\\1', gapSraSummary[i])
        gapExptDescription = "NA"
    else:
        gapExptId = gapExptName = gapExptOrganism = gapExptDescription = "NA"
    if i == 0:
        gapDF = pandas.DataFrame({'Repo':'DbGaP-SRA','Id':gapExptId, 'Name':gapExptName, 'Organism':gapExptOrganism, 'Description':gapExptDescription}, index=[0])
    else:
        gapDF = gapDF.append({'Repo': 'DbGaP-SRA', 'Id': gapExptId, 'Name': gapExptName, 'Organism': gapExptOrganism,
                            'Description': gapExptDescription}, ignore_index=True)
gapDF.to_csv('dbGaPSRAAllTerms.csv', index=False)

###### ENA ######
# Free text search
enaLinkAll = 'https://www.ebi.ac.uk/ena/browser/api/xml/textsearch?domain=sra-study&query=huntington+blood+human'
enaSearchResultAll = getParsedOutput(enaLinkAll, 'xml')
enaSearchCountAll = len(enaSearchResultAll['STUDY_SET']['STUDY']) # 0 results
# No filters available

###### GEO ######
# Uses the Entrez API to serve information like most NCBI resources
# Free text search
geoSearchLinkAll = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=huntington+blood+human'
geoSearchResultAll = getParsedOutput(apiLink=geoSearchLinkAll, respType='xml')
geoSearchCountAll = geoSearchResultAll['eSearchResult']['Count'] # 31 results
# Search with filters, free text for terms without filters
geoSearchLinkFilt = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=huntington+blood+AND+human[orgn]'
geoSearchResultFilt = getParsedOutput(apiLink=geoSearchLinkFilt, respType='xml')
geoSearchCountFilt = geoSearchResultFilt['eSearchResult']['Count'] # 29 results
# The Entrez API retrieves the ids of datasets that match the search and a separate query must be made to get experimental details of those ids
# Get all ids that matched the search
geoSearchAllIdsLink = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&retmax=' + geoSearchCountAll + '&term=huntington+blood+human'
geoSearchAllIdsResult = getParsedOutput(geoSearchAllIdsLink, 'xml')
# Make a loop here with say 50 ids per query
geoIds = geoSearchAllIdsResult['eSearchResult']['IdList']['Id']
geoNumSplits = math.ceil(len(geoIds) / 50)
geoSplitIds = numpy.array_split(geoIds, geoNumSplits)
geoSummaryList = []
for i in range(geoNumSplits):
    tmpIds = geoSplitIds[i]
    tmpGeoSummaryLink = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&id=' + ','.join(tmpIds)
    tmpGeoSummaryResult = getParsedOutput(tmpGeoSummaryLink, 'xml')
    geoSummaryList.append(tmpGeoSummaryResult['eSummaryResult']['DocSum'])
geoSummary = [y for x in geoSummaryList for y in x]
# Experiment details to review
# The XML output from the Entrez API is not parsable so here we use regex to retrieve the relevant details for review
for i in range(len(geoSummary)):
    geoExptId = geoSummary[i]['Id']
    geoExptName = geoSummary[i]['Item'][2]['#text']
    geoExptOrganism = geoSummary[i]['Item'][6]['#text']
    geoExptDescription = geoSummary[i]['Item'][3]['#text']
    if i == 0:
        geoDF = pandas.DataFrame({'Repo':'GEO','Id':geoExptId, 'Name':geoExptName, 'Organism':geoExptOrganism, 'Description':geoExptDescription}, index=[0])
    else:
        geoDF = geoDF.append({'Repo':'GEO', 'Id':geoExptId, 'Name':geoExptName, 'Organism':geoExptOrganism,
                            'Description':geoExptDescription}, ignore_index=True)
geoDF.to_csv('geoAllTerms.csv', index=False)

###### Metabolomics Workbench ######
# Free text search in title
metabWBSearchTitleLink = 'https://www.metabolomicsworkbench.org/rest/study/study_title/huntington/summary'
metabWBSearchTitleResult = getParsedOutput(metabWBSearchTitleLink, 'json') # 1 result

###### OmicsDI ######
# Free text search
omicsdiLinkAll = 'https://www.omicsdi.org:443/ws/dataset/search?query=huntington+blood+human&start=0&size=10&faceCount=0'
omicsdiSearchResultAll = getParsedOutput(apiLink=omicsdiLinkAll, respType='json')
omicsdiCountAll = omicsdiSearchResultAll["count"] # 146 results
# Search with filters, free text for terms without filters
omicsdiLinkAllProper = "https://www.omicsdi.org/ws/dataset/search?query=disease:%22huntington%22%20AND%20TAXONOMY:%229606%22%20AND%20tissue:%22Blood%22"
omicsdiSearchResultAllProper = getParsedOutput(apiLink=omicsdiLinkAllProper, respType='json')
omicsdiCountAllProper = omicsdiSearchResultAllProper["count"] # 0 results
# The API limits the number of search results that can be retrieved but gives an indication of the total number of results
# Here we take the total number of results and split the list into about 100 ids and get information on them iteratively
searchNum = omicsdiCountAll
searchStart = 0
searchSize = 100
searchDone = True
omicsdiAllDatasetsAll = []
while searchDone:
    tmpOmicsdiLink = 'https://www.omicsdi.org:443/ws/dataset/search?query=huntington+blood+human&start=' + str(searchStart) + '&size=' + str(searchSize) +'&faceCount=0'
    tmpOmicsdiSearchResult = getParsedOutput(apiLink=tmpOmicsdiLink, respType='json')
    for i in range(len(tmpOmicsdiSearchResult['datasets'])):
        omicsdiAllDatasetsAll.append(tmpOmicsdiSearchResult['datasets'][i])
    if searchStart + searchSize > searchNum:
        searchDone = False
    else:
        searchStart = searchStart + searchSize
# Experiment details to review
for i in range(len(omicsdiAllDatasetsAll)):
    omicsdiExptId = omicsdiAllDatasetsAll[i]['id']
    omicsdiExptName = omicsdiAllDatasetsAll[i]['title']
    if len(omicsdiAllDatasetsAll[i]['organisms']) != 0:
        omicsdiExptOrganism = omicsdiAllDatasetsAll[i]['organisms'][0]['name']
    omicsdiExptDescription = omicsdiAllDatasetsAll[i]['description']
    if i == 0:
        omicsdiDF = pandas.DataFrame({'Repo':'OmicsDI','Id':omicsdiExptId, 'Name':omicsdiExptName, 'Organism':omicsdiExptOrganism, 'Description':omicsdiExptDescription}, index=[0])
    else:
        omicsdiDF = omicsdiDF.append({'Repo':'OmicsDI', 'Id':omicsdiExptId, 'Name':omicsdiExptName, 'Organism':omicsdiExptOrganism,
                            'Description':omicsdiExptDescription}, ignore_index=True)
omicsdiDF.to_csv('omicsdiAllTerms.csv', index=False)

## PRIDE
# Free text search
# The PRIDE API shows 1000 results per page and the number of results exceeded 1000 so to find out how many results
# there were we query the system interatively until the number of results is less than 1000 which we presume is
# the end of the list
# Final number 5246
prideSearchResultAll = []
for i in range(100):
    prideLinkAll = 'https://www.ebi.ac.uk:443/pride/ws/archive/project/list?query=huntington+blood+human&show=1000&page=' + str(i) + '&order=desc'
    tmpPrideSearchResultAll = getParsedOutput(apiLink=prideLinkAll, respType='json')
    prideSearchCountAll = len(tmpPrideSearchResultAll['list'])
    prideSearchResultAll.append(tmpPrideSearchResultAll['list'])
    print(prideSearchCountAll)
    if prideSearchCountAll == 1000:
        continue
    else:
        break
prideSearchResultAll = [y for x in prideSearchResultAll for y in x]
# Since the number of results is not expected to be over 5k for a rare disease we suspect that the search engine is
# performing an 'OR' search instead of 'AND" with the search terms
# Separate the three queries and do an AND operation
# Disease
prideLinkDisease = 'https://www.ebi.ac.uk:443/pride/ws/archive/project/list?query=huntington&show=100&page=0&order=desc'
prideSearchResultDisease = getParsedOutput(apiLink=prideLinkDisease, respType='json') # 19 results
# only 19 so no need for loop
# Tissue
prideLinkTissue = 'https://www.ebi.ac.uk:443/pride/ws/archive/project/list?query=blood&show=1000&page=0&order=desc'
prideSearchResultTissue = getParsedOutput(apiLink=prideLinkTissue, respType='json')
# exactly 1000 so expanding the query with a for loop
prideSearchResultAllTissue = []
for i in range(100):
    prideLinkAllTissue = 'https://www.ebi.ac.uk:443/pride/ws/archive/project/list?query=blood&show=1000&page=' + str(i) + '&order=desc'
    tmpPrideSearchResultAllTissue = getParsedOutput(apiLink=prideLinkAllTissue, respType='json')
    prideSearchCountAllTissue = len(tmpPrideSearchResultAllTissue['list'])
    prideSearchResultAllTissue.append(tmpPrideSearchResultAllTissue['list'])
    print(prideSearchCountAllTissue)
    if prideSearchCountAllTissue == 1000:
        continue
    else:
        break
prideSearchResultAllTissue = [y for x in prideSearchResultAllTissue for y in x] # 1001 results
# Species
prideLinkSpecies = 'https://www.ebi.ac.uk:443/pride/ws/archive/project/list?query=human&show=100&page=0&order=desc'
# We assume that it is going to be a lot of results because we are searching for human data
prideSearchResultAllSpecies = []
for i in range(100):
    prideLinkAllSpecies = 'https://www.ebi.ac.uk:443/pride/ws/archive/project/list?query=human&show=1000&page=' + str(i) + '&order=desc'
    tmpPrideSearchResultAllSpecies = getParsedOutput(apiLink=prideLinkAllSpecies, respType='json')
    prideSearchCountAllSpecies = len(tmpPrideSearchResultAllSpecies['list'])
    prideSearchResultAllSpecies.append(tmpPrideSearchResultAllSpecies['list'])
    print(prideSearchCountAllSpecies)
    if prideSearchCountAllSpecies == 1000:
        continue
    else:
        break
prideSearchResultAllSpecies = [y for x in prideSearchResultAllSpecies for y in x] # 5019 results
# Do an intersection of all the ids
separateSearchDiseaseIds = []
separateSearchTissueIds = []
separateSearchSpeciesIds = []
for i in range(len(prideSearchResultDisease['list'])):
    separateSearchDiseaseIds.append(prideSearchResultDisease['list'][i]['accession'])
for i in range(len(prideSearchResultAllTissue)):
    separateSearchTissueIds.append(prideSearchResultAllTissue[i]['accession'])
for i in range(len(prideSearchResultAllSpecies)):
    separateSearchSpeciesIds.append(prideSearchResultAllSpecies[i]['accession'])
unqSeparateSearchAllIds = list(set(separateSearchDiseaseIds) & set(separateSearchSpeciesIds)) # 0 results in intersection
# Search with filters
prideLinkAllProper = 'https://www.ebi.ac.uk:443/pride/ws/archive/project/list?show=10&page=0&order=desc&speciesFilter=9606%20&tissueFilter=blood&diseaseFilter=huntington'
prideSearchResultAllProper = getParsedOutput(apiLink=prideLinkAllProper, respType='json')
prideSearchCountAllProper = len(prideSearchResultAllProper['list']) # 0 results

# Save all the API results
sessionDumpFilename = 'globalsaveAllAPIsAllTerms.pkl'
dill.dump_session(sessionDumpFilename)


###### LOAD THE PREVIOUSLY SAVED VARIABLES ######
import dill
dill.load_session('globalsaveAllAPIsAllTerms.pkl')