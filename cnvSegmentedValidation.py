import os
import sys
import requests
import json
import subprocess
import tarfile
import pandas

if ( len(sys.argv) != 3 ):
    print("Usage:\npython3 cnvSegmentedValidation.py [Project Name] [Xena File Path]")
    exit(1)
projectName = sys.argv[1]
# projectName = "CPTAC-3"
xenaFilePath = sys.argv[2]
# xenaFilePath = "/Users/jaimes28/Desktop/gdcData/CPTAC-3/Xena_Matrices/CPTAC-3.segment_cnv_ascat-ngs.tsv"


dataCategory = "copy number variation"
gdcDataType = "Copy Number Segment"
experimentalStrategy = "WGS"
workflowType = "AscatNGS"


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
    for sample in sampleDict:
        # Make success equal to true and make it false if there is an error
        success = True
        # print(f"Sample: {sample}\n\n")

        # Get file id and name so sampleFile name can be constructed
        fileId = sampleDict[sample]["file_id"]
        fileName = sampleDict[sample]["file_name"]
        sampleFile = "gdcFiles/" + fileId + "/" + fileName

        # Create data frame for sample data
        sampleDataDF = pandas.read_csv(sampleFile, sep="\t")

        # Get index of first sample
        firstSampleIndex = xenaDfSampleColumn.index(sample)
        # Find index of last sample and convert it to correct number
        lastSampleIndex = len(xenaDF) - reversedXenaDfSampleColumn.index(sample) - 1

        # Calculate how many rows are in the sample dataframe
        dataCount = lastSampleIndex - firstSampleIndex + 1
        # If data counts don't match then there is missing data so print error
        if dataCount != len(sampleDataDF):
            print(f"Error\nXena matrix does not have all data for sample '{sample}'")
            failedSamples.append(sample)
            continue

        # Loop through each sample
        for sampleIndex in range(firstSampleIndex, lastSampleIndex + 1):
            # Get each row from each dataframe
            xenaRow = xenaDF.iloc[sampleIndex]
            sampleRow = sampleDataDF.iloc[sampleIndex - firstSampleIndex]   # Get proper row in sample dataframe since index is different
            # If any data doesn't match then modify success to show that there was issue and log errors
            if ((sampleRow["Chromosome"] != xenaRow["Chrom"]) or (sampleRow["Start"] != xenaRow["Start"]) or
                    (sampleRow["End"] != xenaRow["End"]) or (sampleRow["Copy_Number"] != xenaRow["value"])):
                success = False
                # Log the errors
                print(f"Error for sample {sample}\n")
                print("Xena Data:\n"
                      f'Chrom: {xenaRow["Chrom"]}\n'
                      f'Start: {xenaRow["Start"]}\n'
                      f'End: {xenaRow["End"]}\n'
                      f'Value: {xenaRow["value"]}\n'
                      )

                print("Sample Data:\n"
                      f'Chrom: {sampleRow["Chromosome"]}\n'
                      f'Start: {sampleRow["Start"]}\n'
                      f'End: {sampleRow["End"]}\n'
                      f'Value: {sampleRow["Copy_Number"]}\n'
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

# print(len(sampleDict), len(xenaSamples))
#
# print(f"Retrieved Samples:\n{sorted(sampleDict)}")
# print(f"Xena Samples:\n{sorted(xenaSamples)}")
if sorted(sampleDict) != sorted(xenaSamples):
    print("ERROR:\nSamples retrieved from GDC do not match those found in xena samples")
    exit(1)

fileIDS = [sampleDict[x]["file_id"] for x in sampleDict]
downloadFiles(fileIDS)

if compare():
    print("Passed")
else:
    print("Fail")
