import os
import sys
import requests
import json
import subprocess
import tarfile
import pandas
import numpy
import difflib
from math import floor, log10

if ( len(sys.argv) != 3 ):
    print("Usage:\npython3 somaticMutationValidation.py [Project Name] [Xena File Path]")
    exit(1)

projectName = sys.argv[1]
# projectName = "HCMI-CMDC"
# xenaFilePath = "/Users/jaimes28/Desktop/gdcData/HCMI-CMDC/Xena_Matrices/HCMI-CMDC.somaticmutation_snv.tsv"
xenaFilePath = sys.argv[2]

dataCategory = "simple nucleotide variation"
gdcDataType = "Masked Somatic Mutation"
experimentalStrategy = "WXS"
workflowType = "Aliquot Ensemble Somatic Variant Merging and Masking"


def round_ForNans(x):
    if( pandas.notna(x) ):
        return numpy.format_float_scientific(x, precision=10)
    else:
        return numpy.nan

# From https://github.com/corriander/python-sigfig/blob/dev/sigfig/sigfig.py
def round_(x, n):
    """Round a float, x, to n significant figures.

	Caution should be applied when performing this operation.
	Significant figures are an implication of precision; arbitrarily
	truncating floats mid-calculation is probably not Good Practice in
	almost all cases.

	Rounding off a float to n s.f. results in a float. Floats are, in
	general, approximations of decimal numbers. The point here is that
	it is very possible to end up with an inexact number:

		roundsf(0.0012395, 3)
		0.00124
	    roundsf(0.0012315, 3)
		0.0012300000000000002

	Basically, rounding in this way probably doesn't do what you want
	it to.
    """
    n = int(n)
    x = float(x)

    if x == 0: return 0

    e = floor(log10(abs(x)) - n + 1)  # exponent, 10 ** e
    shifted_dp = x / (10 ** e)  # decimal place shifted n d.p.
    return round(shifted_dp) * (10 ** e)  # round and revert


def vaf(t_alt_count, t_depth):
    return round_ForNans(t_alt_count / t_depth)


def downloadFiles(fileList):
    jsonPayload = {
        "ids": fileList
    }
    with open("payload.txt", "w") as payloadFile:
        payloadFile.write(str(jsonPayload).replace("\'", "\""))

    print("Downloading from GDC: ")
    subprocess.run(
        ["curl", "-o", "gdcFiles.tar.gz", "--remote-name", "--remote-header-name", "--request", "POST", "--header",
         "Content-Type: application/json", "--data", "@payload.txt", "https://api.gdc.cancer.gov/data"])

    os.system("mkdir -p gdcFiles")
    os.system("tar -xzf gdcFiles.tar.gz -C gdcFiles")


def getXenaSamples(xenaFile):  # get all samples from the xena matrix
    xenaMatrix = pandas.read_csv(xenaFile, sep="\t")
    sampleList = list(xenaMatrix["sample"].unique())
    return sampleList


def getAllSamples(projectName):
    casesEndpt = "https://api.gdc.cancer.gov/cases"
    allSamplesFilter = {
        "op": "and",
        "content": [
            {
                "op": "in",
                "content": {
                    "field": "cases.project.project_id",
                    "value": [
                        projectName
                    ]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "files.analysis.workflow_type",
                    "value": [
                        workflowType
                    ]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "files.data_category",
                    "value": [
                        dataCategory
                    ]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "files.data_type",
                    "value": [
                        gdcDataType
                    ]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "files.experimental_strategy",
                    "value": [
                        experimentalStrategy
                    ]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "files.access",
                    "value": [
                        "open"
                    ]
                }
            }
        ]
    }

    params = {
        "filters": json.dumps(allSamplesFilter),
        "fields": "submitter_sample_ids",
        "format": "json",
        "size": 20000
    }
    response = requests.post(casesEndpt, json=params, headers={"Content-Type": "application/json"})
    responseJson = unpeelJson(response.json())

    allSamples = []
    for caseDict in responseJson:
        for sample in caseDict["submitter_sample_ids"]:
            allSamples.append(sample)
    return allSamples


def unpeelJson(jsonObj):
    jsonObj = jsonObj.get("data").get("hits")
    return jsonObj


def dataTypeSamples(samples):
    filesEndpt = "https://api.gdc.cancer.gov/files"
    # MAKE IT SO THAT FILTER GETS DATA TYPE INSERTED
    dataTypeFilter = {
        "op": "and",
        "content": [
            {
                "op": "in",
                "content": {
                    "field": "cases.project.project_id",
                    "value": [
                        projectName
                    ]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "analysis.workflow_type",
                    "value": [
                        workflowType
                    ]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "data_category",
                    "value": dataCategory
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "data_type",
                    "value": [
                        gdcDataType
                    ]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "experimental_strategy",
                    "value": [
                        experimentalStrategy
                    ]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "cases.samples.submitter_id",
                    "value": samples
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "cases.samples.tissue_type",
                    "value": [
                        "tumor"
                    ]
                }
            }
        ]
    }
    params = {
        "filters": json.dumps(dataTypeFilter),
        "fields": "cases.samples.submitter_id,cases.samples.tissue_type,file_id,file_name",
        "format": "json",
        "size": 20000
    }
    response = requests.post(filesEndpt, json=params, headers={"Content-Type": "application/json"})
    responseJson = unpeelJson(response.json())
    dataTypeDict = {}
    # create seen dict to see how many times a sample has been seen
    seenDict = {}
    seenSamples = []
    for caseDict in responseJson:
        for sample in caseDict["cases"][0]["samples"]:
            sampleName = sample["submitter_id"]
            if sample["tissue_type"] == "Tumor":
                seenDict[sampleName] = seenDict.get(sampleName, 0) + 1
                seenSamples.append(sampleName)
                dataTypeDict[sampleName + "." + str(seenDict[sampleName])] = dict(file_id=caseDict["file_id"],
                                                                                  file_name=caseDict["file_name"])

    return dataTypeDict, list(set(seenSamples))


def xenaDataframe(xenaFile):
    xenaDF = pandas.read_csv(xenaFile, sep="\t")
    for column in xenaDF:
        xenaDF[column] = xenaDF[column].apply(round_ForNans)
    return xenaDF

def nonEmptySamples():
    nonEmpty = []
    allSampleNames = []
    for sample in sampleDict:
        fileId = sampleDict[sample]["file_id"]
        fileName = sampleDict[sample]["file_name"]
        sampleFile = "gdcFiles/" + fileId + "/" + fileName
        normalSampleName = sample[:sample.index(".")]
        sampleDataDF = pandas.read_csv(sampleFile, sep="\t", skiprows=7)
        if (len(sampleDataDF.index) == 0):
            continue
        if (normalSampleName not in nonEmpty):
            nonEmpty.append(normalSampleName)
            allSampleNames.append(sample)

    return nonEmpty, allSampleNames


def sampleDataframe():
    dataFrame = pandas.DataFrame()
    # Create a dataframe for all the samples retrieved
    for sample in sampleNames:
        fileId = sampleDict[sample]["file_id"]
        fileName = sampleDict[sample]["file_name"]
        sampleFile = "gdcFiles/" + fileId + "/" + fileName
        normalSampleName = sample[:sample.index(".")]
        # Create data frame for sample data
        sampleDataDF = pandas.read_csv(sampleFile, sep="\t", skiprows=7)
        sampleDataDF.rename(columns={'Hugo_Symbol': 'gene'}, inplace=True)
        sampleDataDF.rename(columns={'Chromosome': 'chrom'}, inplace=True)
        sampleDataDF.rename(columns={'Start_Position': 'start'}, inplace=True)
        sampleDataDF.rename(columns={'End_Position': 'end'}, inplace=True)
        sampleDataDF.rename(columns={'Reference_Allele': 'ref'}, inplace=True)
        sampleDataDF.rename(columns={'Tumor_Seq_Allele2': 'alt'}, inplace=True)
        sampleDataDF.rename(columns={'HGVSp_Short': 'Amino_Acid_Change'}, inplace=True)
        sampleDataDF.rename(columns={'Consequence': 'effect'}, inplace=True)

        sampleDataDF["dna_vaf"] = sampleDataDF.apply(lambda x: vaf(x.t_alt_count, x.t_depth), axis=1)

        sampleDataDF = sampleDataDF.loc[:,
                       ["gene", "chrom", "start", "end", "ref", "alt", "Tumor_Sample_Barcode", "Amino_Acid_Change",
                        "effect",
                        "callers", "dna_vaf"]]
        sampleDataDF = sampleDataDF[
            ["gene", "chrom", "start", "end", "ref", "alt", "Tumor_Sample_Barcode", "Amino_Acid_Change", "effect",
             "callers", "dna_vaf"]]
        sampleDataDF.insert(0, "sample", normalSampleName)
        sampleDataDF["dna_vaf"] = sampleDataDF["dna_vaf"].apply(round_ForNans)

        dataFrame = pandas.concat([dataFrame, sampleDataDF])
    return dataFrame


xenaSamples = getXenaSamples(xenaFilePath)

allSamples = getAllSamples(projectName)
sampleDict, seenSamples = dataTypeSamples(allSamples)
xenaDF = xenaDataframe(xenaFilePath)

fileIDS = [sampleDict[x]["file_id"] for x in sampleDict]
downloadFiles(fileIDS)

validSamples, sampleNames = nonEmptySamples()
validSamples = list(set(validSamples))

if sorted(validSamples) != sorted(xenaSamples):
    print("Samples retrieved from GDC and not in Xena Dataframe")
    print([x for x in seenSamples if x not in xenaSamples])
    print("Samples retrieved from Xena Dataframe and not in GDC retrieved samples")
    print([x for x in xenaSamples if x not in seenSamples])
    exit(1)

# create dataframe for samples
sampleDf = sampleDataframe()
xenaDF.sort_values(by=sorted(xenaDF), inplace=True)
# sort sample dataframe as well
sampleDf.sort_values(by=sorted(xenaDF), inplace=True)

# then reset index ordering for each one
xenaDF.reset_index(inplace=True, drop=True)
sampleDf.reset_index(inplace=True, drop=True)

with open("sampleDF.csv", "w") as sampleFile:
    sampleDf.to_csv(sampleFile)

with open("xenaDF.csv", "w") as xenaDfFile:
    xenaDF.to_csv(xenaDfFile)


try:
    pandas.testing.assert_frame_equal(sampleDf, xenaDF, check_dtype=False)
    print("Passed")
except AssertionError:
    with open("sampleDF.csv", "r") as sampleFile:
        with open("xenaDF.csv", "r") as xenaDfFile:
            # if they are not equal then output diff of both files
            with open("diff.txt", "w") as diffFile:
                diffFile.writelines(difflib.unified_diff(sampleFile.readlines(), xenaDfFile.readlines(),
                                                         fromfile="sampleDF.csv", tofile="xenaDF.csv"))
