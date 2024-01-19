import os
import sys
import requests
import json
import subprocess
import tarfile
import pandas
import numpy

if ( len(sys.argv) != 3 ):
    print("Usage:\npython3 somaticMutationValidation.py [Project Name] [Xena File Path]")
    exit(1)

projectName = sys.argv[1]
# projectName = "HCMI-CMDC"
# xenaFilePath = "/Users/jaimes28/Downloads/HCMI-CMDC.somaticmutation_snv (1).tsv"
xenaFilePath = sys.argv[2]

dataCategory = "simple nucleotide variation"
gdcDataType = "Masked Somatic Mutation"
experimentalStrategy = "WXS"
workflowType = "Aliquot Ensemble Somatic Variant Merging and Masking"



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

    gdcTarfile = tarfile.open("gdcFiles.tar.gz")
    gdcTarfile.extractall("gdcFiles")
    gdcTarfile.close()


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
            },
            {
                "op": "in",
                "content": {
                    "field": "access",
                    "value": [
                        "open"
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
    for caseDict in responseJson:
        for sample in caseDict["cases"][0]["samples"]:
            if sample["tissue_type"] == "Tumor":
                dataTypeDict[sample["submitter_id"]] = dict(file_id=caseDict["file_id"],
                                                            file_name=caseDict["file_name"])

    return dataTypeDict


def xenaDataframe(xenaFile):
    xenaDF = pandas.read_csv(xenaFile, sep="\t")
    xenaDF = xenaDF.round(10)
    return xenaDF


def compare():
    # Get a list of xena data frame column to see where sample name first occurs
    xenaDfSampleColumn = xenaDF["sample"].tolist()
    # get reversed list to find last occurence of sample in the column
    reversedXenaDfSampleColumn = list(reversed(xenaDfSampleColumn))

    # Create empty list to collect all the samples that had errors
    failedSamples = []

    # Go through each sample
    for sample in xenaSamples:
        # Make success equal to true and make it false if there is an error
        success = True
        # print(f"Sample: {sample}\n\n")

        # Get file id and name so sampleFile name can be constructed
        try:
            fileId = sampleDict[sample]["file_id"]
        except KeyError:
            print(sample)
            continue
        fileName = sampleDict[sample]["file_name"]
        sampleFile = "gdcFiles/" + fileId + "/" + fileName

        # Create data frame for sample data
        sampleDataDF = pandas.read_csv(sampleFile, sep="\t", skiprows=7)

        # Get index of first sample
        firstSampleIndex = xenaDfSampleColumn.index(sample)
        # Find index of last sample and convert it to correct number
        lastSampleIndex = len(xenaDF) - reversedXenaDfSampleColumn.index(sample) - 1

        # Calculate how many rows are in the sample dataframe
        dataCount = lastSampleIndex - firstSampleIndex + 1
        # If data counts don't match then there is missing data so print error
        if dataCount != len(sampleDataDF):
            print(f"Error\nXena matrix does not have all data for sample '{sample}'")
            continue

        # Loop through each sample
        for sampleIndex in range(firstSampleIndex, lastSampleIndex + 1):
            # Get each row from each dataframe
            xenaRow = xenaDF.iloc[sampleIndex]
            sampleRow = sampleDataDF.iloc[sampleIndex - firstSampleIndex]   # Get proper row in sample dataframe since index is different
            dna_vaf = numpy.round(sampleRow["t_alt_count"] / sampleRow["t_depth"], 10)
            # If any data doesn't match then modify success to show that there was issue and log errors
            if ((sampleRow["Hugo_Symbol"] != xenaRow["gene"] and (not pandas.isna(sampleRow["Hugo_Symbol"]) and not pandas.isna(xenaRow["gene"]))) or
                    (sampleRow["Chromosome"] != xenaRow["chrom"] and (not pandas.isna(sampleRow["Chromosome"]) and not pandas.isna(xenaRow["chrom"]))) or
                    (sampleRow["Start_Position"] != xenaRow["start"] and (not pandas.isna(sampleRow["Start_Position"]) and not pandas.isna(xenaRow["start"]))) or
                    (sampleRow["End_Position"] != xenaRow["end"] and (not pandas.isna(sampleRow["End_Position"]) and not pandas.isna(xenaRow["end"]))) or
                    (sampleRow["Reference_Allele"] != xenaRow["ref"] and (not pandas.isna(sampleRow["Reference_Allele"]) and not pandas.isna(xenaRow["ref"]))) or
                    (sampleRow["Tumor_Seq_Allele2"] != xenaRow["alt"] and (not pandas.isna(sampleRow["Tumor_Seq_Allele2"]) and not pandas.isna(xenaRow["alt"]))) or
                    (sampleRow["Tumor_Sample_Barcode"] != xenaRow["Tumor_Sample_Barcode"] and (not pandas.isna(sampleRow["Tumor_Sample_Barcode"]) and not pandas.isna(xenaRow["Tumor_Sample_Barcode"]))) or
                    (sampleRow["HGVSp_Short"] != xenaRow["Amino_Acid_Change"] and (not pandas.isna(sampleRow["HGVSp_Short"]) and not pandas.isna(xenaRow["Amino_Acid_Change"]))) or
                    (sampleRow["Consequence"] != xenaRow["effect"] and (not pandas.isna(sampleRow["Consequence"]) and not pandas.isna(xenaRow["effect"]))) or
                    (sampleRow["callers"] != xenaRow["callers"] and (not pandas.isna(sampleRow["callers"]) and not pandas.isna(xenaRow["callers"]))) or
                    (dna_vaf != xenaRow["dna_vaf"])):
                success = False
                # Log the errors
                print(f"Error for sample {sample}\n")
                print("Xena Data:\n"
                      f'Gene: {xenaRow["gene"]}\n'
                      f'Chrom: {xenaRow["chrom"]}\n'
                      f'Start: {xenaRow["start"]}\n'
                      f'End: {xenaRow["end"]}\n'
                      f'Ref: {xenaRow["ref"]}\n'
                      f'Alt: {xenaRow["alt"]}\n'
                      f'Tumor Sample Barcode: {xenaRow["Tumor_Sample_Barcode"]}\n'
                      f'Amino Acid Change: {xenaRow["Amino_Acid_Change"]}\n'
                      f'Effect: {xenaRow["effect"]}\n'
                      f'Callers: {xenaRow["callers"]}\n'
                      f'Dna Vaf: {xenaRow["dna_vaf"]}\n'
                      )

                print("Sample Data:\n"
                      f'Gene: {sampleRow["Hugo_Symbol"]}\n'
                      f'Chrom: {sampleRow["Chromosome"]}\n'
                      f'Start: {sampleRow["Start_Position"]}\n'
                      f'End: {sampleRow["End_Position"]}\n'
                      f'Ref: {sampleRow["Reference_Allele"]}\n'
                      f'Alt: {sampleRow["Tumor_Seq_Allele2"]}\n'
                      f'Tumor Sample Barcode: {sampleRow["Tumor_Sample_Barcode"]}\n'
                      f'Amino Acid Change: {sampleRow["HGVSp_Short"]}\n'
                      f'Effect: {sampleRow["Consequence"]}\n'
                      f'Callers: {sampleRow["callers"]}\n'
                      f'Dna Vaf: {dna_vaf}\n'
                      )
        # If there was an error then add it to list of failed samples
        if not success:
            failedSamples.append(sample)
    # If there are no failed samples then return True
    if len(failedSamples) == 0:
        return True
    # Otherwise print all the failed samples and return False
    else:
        print(f"Failed Samples:\n {failedSamples}")
        return False



xenaSamples = getXenaSamples(xenaFilePath)

allSamples = getAllSamples(projectName)
sampleDict = dataTypeSamples(allSamples)
xenaDF = xenaDataframe(xenaFilePath)

modifiedSamples = [sample for sample in sampleDict if sample in xenaSamples]
# print(len(modifiedSamples), len(xenaSamples))


# print(f"Retrieved Samples:\n{sorted(sampleDict)}")
# print(f"Xena Samples:\n{sorted(xenaSamples)}")

# UNCOMMENT AFTER
# if sorted(sampleDict) != sorted(xenaSamples):
#     print("ERROR:\nSamples retrieved from GDC do not match those found in xena samples\nMissing:")
#     print([sample for sample in sampleDict if (sample not in xenaSamples)])
#    exit(1)

fileIDS = [sampleDict[x]["file_id"] for x in sampleDict]
downloadFiles(fileIDS)

if compare():
    print("Passed")
else:
    print("Fail")
