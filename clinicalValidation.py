import sys
import requests
import json
import pandas
import numpy
import pandas
import ast

GDC_DROPPED_FIELDS = {
    'cases.aliquot_ids',
    'cases.submitter_aliquot_ids',
    'cases.created_datetime',
    'cases.sample_ids',
    'diagnosis_ids',    # cases.diagnoses.diagnosis_id
    'cases.submitter_sample_ids',
    'submitter_diagnosis_ids', # no submitter diagnosis id
    'cases.updated_datetime',
    'cases.index_date', #does not exist
    'cases.state',
    'cases.project.dbgap_accession_number',
    'cases.project.releasable', # no project.releasable
    'cases.project.state',
    'cases.project.program.dbgap_accession_number',
    'cases.project.program.program_id',
    'cases.project.released',
    'cases.diagnoses.created_datetime',
    'cases.diagnoses.updated_datetime',
    'cases.diagnoses.state',
    'cases.diagnoses.submitter_id',
    'cases.diagnoses.diagnosis_id',
    'cases.demographic.submitter_id',
    'cases.demographic.created_datetime',
    'cases.demographic.demographic_id',
    'cases.demographic.updated_datetime',
    'cases.demographic.state',
    'cases.submitter_slide_ids',
    'cases.submitter_analyte_ids',
    'cases.follow_ups', #Doesn't exist
    'cases.portion_ids',
    'cases.submitter_portion_ids',
    'case.slide_ids',
    'cases.analyte_ids',
    'diagnoses',  # Doesn't exist  diagnose._____
    'diagnoses.treatments', # Doesn't exist diagnoses.treatments._______
    'cases.family_histories.updated_datetime',
    'cases.family_histories.submitter_id',
    'cases.family_histories.state',
    'cases.family_histories.created_datetime',
    'cases.family_histories.family_history_id',
    'cases.exposures.submitter_id',
    'cases.exposures.created_datetime',
    'cases.exposures.updated_datetime',
    'cases.exposures.exposure_id',
    'cases.exposures.state',
    'cases.samples.created_datetime',
    'cases.samples.updated_datetime',
    'cases.samples.state',
    'samples.portions',  #cases.samples.portions_______
}
casesEndpt = "https://api.gdc.cancer.gov/cases"

if (len(sys.argv) != 3):
    print("Usage:\npython3 clinicalValidation.py [Project Name] [Xena File Path]")
    exit(1)

projectName = sys.argv[1]
# projectName = "TARGET-AML"
# xenaFilePath = "/Users/jaimes28/Desktop/gdcData/TARGET-AML/Xena_Matrices/TARGET-AML.clinical.tsv"
xenaFilePath = sys.argv[2]


def unpack_dict(d, parent_key='', sep='.'):
    """
    Unpacks a dictionary recursively, handling lists of dictionaries,
    and returns a new dictionary with composite keys and lists of values.
    Leaves specified keys unchanged and handles single-item lists by
    converting them to the item itself.

    :param d: The dictionary to unpack.
    :param parent_key: The base key string used to build composite keys.
    :param sep: Separator used in composite keys.
    :param special_keys: Set of keys to be left unchanged.
    """

    items = {}

    # Iterate through the dictionary
    for k, v in d.items():
        new_key = f"{parent_key}{sep}{k}" if parent_key else k

        if k == "samples":
            items[new_key] = v
        elif isinstance(v, dict):
            # Recurse into nested dictionaries
            items.update(unpack_dict(v, new_key, sep))
        elif isinstance(v, list):
            if len(v) == 1:
                # If the list contains a single item, use the item itself
                v = v[0]
                if isinstance(v, dict):
                    # If it's a single dictionary, unpack its values directly
                    items.update(unpack_dict(v, new_key, sep))
                else:
                    # If it's a single non-dict item, add it directly
                    items[new_key] = v
            elif all(isinstance(i, dict) for i in v):
                # unpack this before going into it?
                v = [unpack_dict(x) for x in v]
                # Handle lists of dictionaries
                all_keys = set()
                # maybe don't include keys that have list of dicts as values?
                for item in v:
                    all_keys.update(item.keys())

                # Initialize new dictionary with composite keys and empty lists
                for sub_key in all_keys:
                    composite_key = f"{new_key}{sep}{sub_key}"
                    items[composite_key] = []

                # Populate the lists with values from each dictionary in the list
                for sub_key in all_keys:
                    composite_key = f"{new_key}{sep}{sub_key}"
                    for item in v:
                        retrievedValue = item.get(sub_key, '')
                        items[composite_key].append(retrievedValue if retrievedValue is not None else '')

                    # if sub_key.startswith(("treatments", "pathology_details", "annotations")) and any(isinstance(i, list) for i in items[composite_key]):
                    #     organizedList = []
                    #     for entry in items[composite_key]:
                    #         if isinstance(entry, list):
                    #             organizedList.extend(entry)
                    #     items[composite_key] = organizedList
            else:
                # Handle non-dict elements in lists (if necessary)
                items[new_key] = v
        else:
            # Directly add the value if it's not a nested dictionary or list of dictionaries
            items[new_key] = v

    return items


def customMin(minList):
    convertedList = []
    for value in minList:
        try:
            convertedList.append(float(value))
        except:
            convertedList.append(float("inf"))
    return min(convertedList)


def general_compare(val1, val2):
    # Normalize "None" and "NaN" string representations
    if isinstance(val1, str):
        val1_lower = val1.lower()
        if val1_lower == "none":
            val1 = None
        elif val1_lower == "nan":
            val1 = float('nan')
    if isinstance(val2, str):
        val2_lower = val2.lower()
        if val2_lower == "none":
            val2 = None
        elif val2_lower == "nan":
            val2 = float('nan')

    # Handle None and NaN equivalence
    if (val1 is None and isinstance(val2, float) and numpy.isnan(val2)) or \
            (val2 is None and isinstance(val1, float) and numpy.isnan(val1)):
        return True

    # Handle None values
    if val1 is None or val2 is None:
        return val1 is val2

    # Handle NaN values
    if isinstance(val1, float) and isinstance(val2, float):
        if numpy.isnan(val1) and numpy.isnan(val2):
            return True

    # Convert string representation of a list to an actual list if necessary
    if isinstance(val1, str) and val1.startswith('[') and val1.endswith(']'):
        try:
            val1 = ast.literal_eval(val1)
        except (ValueError, SyntaxError):
            return False
    if isinstance(val2, str) and val2.startswith('[') and val2.endswith(']'):
        try:
            val2 = ast.literal_eval(val2)
        except (ValueError, SyntaxError):
            return False

    # Ensure both values are lists before comparing
    if isinstance(val1, list) and isinstance(val2, list):
        if len(val1) != len(val2):
            return False
        for v1, v2 in zip(val1, val2):
            if not general_compare(v1, v2):
                return False
        return True

    # Attempt to convert both values to floats for comparison
    try:
        return numpy.format_float_scientific(float(val1), precision=8) == numpy.format_float_scientific(float(val2),
                                                                                                  precision=8)
    except (ValueError, TypeError):
        # If conversion to float fails, compare values as strings
        return str(val1).lower() == str(val2).lower()


def flatten_json(json_obj, prefix=""):
    flattened_json = {}
    for key, value in json_obj.items():
        if isinstance(value, dict):
            flattened_json.update(flatten_json(value, prefix + key + "."))
        else:
            flattened_json[prefix + key] = value
    return flattened_json


def update_dict_without_overwriting(original_dict, new_dict):
    """
    Update the original dictionary with the new dictionary’s values,
    appending to the existing values without overwriting.
    :param original_dict: The dictionary to update.
    :param new_dict: The dictionary with new values.
    :return: The updated dictionary.
    """
    for key, value in new_dict.items():
        if key in original_dict:
            if value == original_dict[key]:
                continue
            # If the key exists and both values are lists, extend the list
            if isinstance(original_dict[key], list) and isinstance(value, list):
                original_dict[key].extend(value)
            # If the key exists and both values are sets, update the set
            elif isinstance(original_dict[key], set) and isinstance(value, set):
                original_dict[key].update(value)
            # If the key exists and both values are strings, concatenate the strings
            elif isinstance(original_dict[key], str) and isinstance(value, str):
                original_dict[key] = [original_dict[key]]
                original_dict[key].append(value)
            # If the key exists and both values are dictionaries, update the dictionary recursively
            elif isinstance(original_dict[key], dict) and isinstance(value, dict):
                update_dict_without_overwriting(original_dict[key], value)
            # Otherwise, convert both values to a list and append the new value
            else:
                if not isinstance(original_dict[key], list):
                    original_dict[key] = [original_dict[key]]
                original_dict[key].append(value)
        else:
            original_dict[key] = value
    return original_dict


def getAllSamples(projectName):
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
            }
        ]
    }

    params = {
        "filters": json.dumps(allSamplesFilter),
        "fields": "submitter_sample_ids",
        "format": "json",
        "size": 2000000
    }
    response = requests.post(casesEndpt, json=params, headers={"Content-Type": "application/json"})
    responseJson = unpeelJson(response.json())

    allSamples = []
    for caseDict in responseJson:
        for sample in caseDict["submitter_sample_ids"]:
            allSamples.append(sample)
    return allSamples


def availableFields():
    response = requests.get("https://api.gdc.cancer.gov/files/_mapping")
    availableFields = response.json()["fields"]
    wantedFields = []
    for field in availableFields:
        if field not in GDC_DROPPED_FIELDS and not field.startswith('cases.summary.') and field.startswith("cases.") and not field.startswith("cases.follow_ups."):
            wantedFields.append(field)
    return wantedFields


def unpeelJson(jsonObj):
    jsonObj = jsonObj.get("data").get("hits")
    return jsonObj


def getFieldData(fields):
    projectFields = [x.removeprefix("cases.project.") for x in fields if "cases.project." in x]
    for x in fields:
        if "cases.project." in x:
            fields.remove(x)

    sampleFields = [x.removeprefix("cases.") for x in fields if "cases.samples." in x]
    for x in fields:
        if "cases.samples." in x:
            fields.remove(x)

    projectFilter = {
        "op": "and",
        "content": [
            {
                "op": "in",
                "content": {
                    "field": "project_id",
                    "value": [
                        projectName
                    ]
                }
            }
        ]
    }
    caseFilter = {
        "op": "and",
        "content": [
            {
                "op": "in",
                "content": {
                    "field": "project.project_id",
                    "value": [
                        projectName
                    ]
                }
            }
        ]
    }
    fileFilter = {
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
            }
        ]
    }

    caseParams = {
        "filters": json.dumps(caseFilter),
        "fields": ",".join(sampleFields) + ",case_id",
        "format": "json",
        "size": 2000000
    }

    response = requests.post("https://api.gdc.cancer.gov/cases/", json=caseParams, headers={"Content-Type": "application/json"})
    sampleResponseJson = unpeelJson(response.json())
    sampleResponseDictionary = {case["case_id"]: case for case in sampleResponseJson}
    for case in sampleResponseDictionary:
        for key in list(sampleResponseDictionary[case].keys()):
            if key != "samples":
                del sampleResponseDictionary[case][key]

    projectParams = {
        "filters": json.dumps(projectFilter),
        "fields": ",".join(projectFields),
        "format": "json",
        "size": 2000000
    }

    response = requests.post("https://api.gdc.cancer.gov/projects/", json=projectParams, headers={"Content-Type": "application/json"})
    projectResponseJson = unpeelJson(response.json())
    projectData = {"project": projectResponseJson[0]}

    callOneFields = fields[(len(fields)//2):]

    params = {
        "filters": json.dumps(fileFilter),
        "fields": ",".join(callOneFields) + ",cases.case_id" + ',cases.diagnoses.treatments.treatment_id'+ ",cases.diagnoses.diagnosis_id",
        "format": "json",
        "size": 2000000
    }

    response = requests.post("https://api.gdc.cancer.gov/files/", json=params, headers={"Content-Type": "application/json"})
    response1Json = unpeelJson(response.json())
    response1Dictionary = {case["cases"][0]["case_id"]: case["cases"][0] for case in response1Json}

    callTwoFields = fields[:(len(fields)//2)]

    params = {
        "filters": json.dumps(fileFilter),
        "fields": ",".join(callTwoFields) + ",cases.case_id" + ',cases.diagnoses.treatments.treatment_id' + ",cases.diagnoses.diagnosis_id",
        "format": "json",
        "size": 2000000
    }
    response = requests.post("https://api.gdc.cancer.gov/files/", json=params, headers={"Content-Type": "application/json"})
    response2Json = unpeelJson(response.json())
    response2Dictionary = {case["cases"][0]["case_id"]: case["cases"][0] for case in response2Json}

    combinedData = update_dict_without_overwriting(response1Dictionary, response2Dictionary)
    del projectData["project"]["id"]

    for caseId in combinedData:
        combinedData[caseId].update(projectData)
        combinedData[caseId].update(sampleResponseDictionary[caseId])
    return combinedData


# ***
def validSamples(caseData):
    fileFilter = {
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
                    "field": "access",
                    "value": [
                        'open'
                    ]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "files.data_category",
                    "value": [
                        'transcriptome profiling', 'proteome profiling', 'dna methylation', 'copy number variation', 'simple nucleotide variation'
                    ]
                }
            }
        ]
    }
    params = {
        "filters": json.dumps(fileFilter),
        "fields": 'data_category,cases.samples.submitter_id,cases.samples.tissue_type',
        "format": "json",
        "size": 2000000
    }

    response = requests.post("https://api.gdc.cancer.gov/files/", json=params, headers={"Content-Type": "application/json"})
    response = unpeelJson(response.json())
    keepSamples = []
    totalSamples = []
    for case in caseData:
        case = caseData[case]
        for sample in case["samples"]:
            totalSamples.append(sample["submitter_id"])



    for file in response:
        samples = file["cases"][0]["samples"]
        # Get samples associated with transcriptome, proteome, and/or methylation file(s)
        if file['data_category'] in ['Transcriptome Profiling', 'Proteome Profiling', 'DNA Methylation']:
            for sample in samples:
                submitter_id = sample['submitter_id']
                if submitter_id not in keepSamples:
                    keepSamples.append(submitter_id)
        # Get 'Tumor' samples associated with copy number variation and/or simple nucleotide variation file(s)
        if file['data_category'] in ['Copy Number Variation', 'Simple Nucleotide Variation']:
            for sample in samples:
                submitter_id = sample['submitter_id']
                if sample['tissue_type'] == 'Tumor' and submitter_id not in keepSamples:
                    keepSamples.append(submitter_id)

    rightSamples = 0
    for caseID, caseInfo in caseData.items():
        sampleIndex = 0
        while sampleIndex < (len(caseInfo["samples"])):
            sample = caseInfo["samples"][sampleIndex]
            if sample["submitter_id"] not in keepSamples:
                caseInfo["samples"].pop(sampleIndex)
                sampleIndex -= 1
            else:
                rightSamples += 1
            sampleIndex += 1
    x=5


# *****
def formatDiagnosis(caseData):
    for case in caseData:
        case = caseData[case]
        if "diagnoses" in case:
            diagnosisTempDict = {}
            for diagnosis in case["diagnoses"]:
                id = diagnosis["diagnosis_id"]
                if id not in diagnosisTempDict:
                    del diagnosis["diagnosis_id"]
                    diagnosisTempDict[id] = diagnosis
                else:
                    del diagnosis["diagnosis_id"]
                    update_dict_without_overwriting(diagnosisTempDict[id], diagnosis)
                # diagnosisTempDict = update_dict_without_overwriting(diagnosisTempDict, diagnosis)
            case["diagnoses"] = [diagnosisTempDict[diagnosis] for diagnosis in diagnosisTempDict]


# ***
def formatTreatments(caseData):
    for case in caseData:
        case = caseData[case]
        if "diagnoses" in case:
            for diagnosis in case["diagnoses"]:

                if "treatments" in diagnosis:
                    treatmentTempDict = {}
                    for treatment in diagnosis["treatments"]:
                        if treatment["treatment_id"] not in treatmentTempDict:
                            treatmentTempDict[treatment["treatment_id"]] = treatment
                        else:
                            update_dict_without_overwriting(treatmentTempDict[treatment["treatment_id"]], treatment )

                    diagnosis["treatments"] = [treatmentTempDict[treatment] for treatment in treatmentTempDict]


def formatPathology(caseData):
    for case in caseData:
        case = caseData[case]
        pathologyKeys = set()
        if "diagnoses" in case:
            for diagnosis in case["diagnoses"]:
                if "pathology_details" in diagnosis:
                    tempPathologyDict = {x["pathology_detail_id"]: x for x in case["diagnoses"]["pathology_details"]}
                    case["diagnoses"]["pathology_details"] = tempPathologyDict
                    for pathology in case["diagnoses"]["pathology_details"]:
                        pathology = case["diagnoses"]["pathology_details"][pathology]
                        pathologyKeys.update(set(pathology.keys()))

            pathologyDict = {key: [] for key in pathologyKeys}

            for diagnosis in case["diagnoses"]:
                if "pathology_details" in diagnosis:
                    for pathology in case["diagnoses"]["pathology_details"]:
                        pathology = case["diagnoses"]["pathology_details"][pathology]
                        for key in pathologyDict:
                            pathologyDict[key].append(pathology.get(key, ''))

    # x = [sample["submitter_id"] for case in caseData for sample in caseData[case]["samples"]]
    # df = pandas.read_csv("/Users/jaimes28/Desktop/gdcData/CGCI-HTMCP-DLBCL/Xena_Matrices/CGCI-HTMCP-DLBCL.clinical.tsv", sep="\t")
    # samples = list(df["sample"])
    # t = 5


def formatAnnotations(caseData):
    for case in caseData:
        case = caseData[case]
        annotationKeys = set()
        if "annotations" in case:
            tempAnnotations = {x["annotation_id"]: x for x in case["annotations"]}
            case["annotations"] = tempAnnotations

            for annotation in case["annotations"]:
                annotation = case["annotations"][annotation]
                annotationKeys.update(annotation.keys())

            annotationDict = {key: [] for key in annotationKeys}

            for annotation in case["annotations"]:
                annotation = case["annotations"][annotation]
                for key in annotation:
                    annotationDict[key].append(annotation.get(key, ''))
            if len(case["annotations"]) == 1:
                annotationDict = {key: annotationDict[key][0] for key in annotationDict}


# ***
# def unpack(caseData):
#     for caseName in caseData:
#         # Flatten case
#         case = flatten_json(caseData[caseName])
#         if isinstance(case["diagnoses"], list) and len(case["diagnoses"]) == 1:
#             case["diagnoses"] = case["diagnoses"][0]
#             case = flatten_json(case)
#             x=5
#         elif isinstance(case["diagnoses"], list) and len(case["diagnoses"]) > 1:
#             allKeys = {key for dictionary in case["diagnoses"] for key in dictionary if key != "pathology_details" and key != "treatments"}
#             tempDict = {"diagnoses." + key: [] for key in allKeys}
#             for dictionary in case["diagnoses"]:
#                 for key in allKeys:
#                     tempDict["diagnoses." + key].append(dictionary.get(key, '') if dictionary.get(key, '') is not None else '')
#             pathologyHolder = [x for x in case["diagnoses"] if "pathology_details" in x]
#             treatmentHolder = [x for x in case["diagnoses"] if "treatments" in x]
#             tempDict.update({"diagnoses.treatments": treatmentHolder[0]["treatments"]})
#             tempDict.update({"diagnoses.pathology_details": pathologyHolder[0]["pathology_details"]})
#             case.update(tempDict)
#             del case["diagnoses"]
#
#         if "diagnoses.treatments" in case:
#             treatmentTempDict = {}
#             for treatment in case["diagnoses.treatments"]:
#                 if treatment["treatment_id"] not in treatmentTempDict:
#                     treatmentTempDict[treatment["treatment_id"]] = treatment
#                 else:
#                     update_dict_without_overwriting(treatmentTempDict[treatment["treatment_id"]], treatment)
#
#             case["diagnoses.treatments"] = [treatmentTempDict[treatment] for treatment in treatmentTempDict]
#
#
#
#
#         # collect all keys to remove from caseData
#         keysToDelete = []
#         # create new dictionary to collect unpacked data
#         unpackedDict = {}
#
#         for key, value in case.items():
#             # if we have a list of dictionaries (and not samples)
#             if value is None:
#                 keysToDelete.append(key)
#                 continue
#             if isinstance(value, list) and len(value) == 1 and not isinstance(value[0], dict):
#                 case[key] = value[0]
#                 continue
#             if key != "samples" and isinstance(value, list) and len(value) > 0 and isinstance(value[0], dict):
#                 # Add to the list the unpacked key that we should delete from original caseData dict
#                 keysToDelete.append(key)
#                 # Create temp value list to hold value with Null values removed
#                 tempValue = []
#                 # add each dict without null values to list
#                 for dictionary in value:
#                     tempDict = {key: value for key, value in dictionary.items() if value is not None}
#                     tempValue.append(tempDict)
#                 # assign value to this
#                 value = tempValue
#
#                 # if value only has one dict then no need to format it like a list
#                 if len(value) == 1:
#                     # iterate through each entry in the dict and add it to the unpacked dict and add parent dict at
#                     # beginning of key
#                     for entry in value[0]:
#                         unpackedDict[key + "." + entry] = value[0][entry]
#                     # go onto next key value pair
#                     continue
#
#                 # if there are multiple dicts collect all keys across each dict in a set
#                 allKeys = {key for dictionary in value for key in dictionary}
#                 # create a temp dict to hold all these keys with an empty list as the value
#                 tempDict = {key + "." + nestedKey: [] for nestedKey in allKeys}
#                 # for each dictionary in value
#                 for dictionary in value:
#                     # get every key from the dictionary and append it to the temp dict
#                     for key2 in allKeys:
#                         # use .get to handle case where key is not in that dictionary (Key Error) and return an empty
#                         # string to append
#                         tempDict[key + "." + key2].append(dictionary.get(key2, ''))
#                 unpackedDict.update(tempDict)
#         case.update(unpackedDict)
#         for key in keysToDelete:
#             del case[key]
#         caseData[caseName] = case
#     print("here")


def processSamples(caseData):
    sampleOrientedData = {}
    for caseName in caseData:
        case = caseData[caseName]
        samples = case["samples"]
        for sample in samples:
            sampleDict = {"samples." + key: value for key, value in sample.items() if value is not None}
            sampleDict.update(case)
            del sampleDict["samples"]
            sampleOrientedData[sampleDict["samples.submitter_id"]] = sampleDict

    return sampleOrientedData


def compare():
    for i in range(0, len(myDF)):
        myDFRow = myDF.iloc[i]
        dataFrameRow = dataFrame.iloc[i]
        for column in myDF:
            if column.endswith((".treatments.diagnoses", ".annotations.diagnoses", ".pathology_details.diagnoses")):
                continue
            xenaCell = dataFrameRow[column]
            myCell = myDFRow[column]
            if pandas.isna(xenaCell) and pandas.isna(myCell):
                continue
            if not general_compare(xenaCell, myCell):
                print("Error")
                print(f"Xena Value: {xenaCell}")
                print(f"Retrieved Value: {myCell}")
                exit(1)
    print("Passed")


samples = getAllSamples(projectName)
wantedFields = availableFields()
caseData = getFieldData(wantedFields)
formatDiagnosis(caseData)
formatTreatments(caseData)

for id in caseData:
    caseData[id] = unpack_dict(caseData[id])

validSamples(caseData)
sampleOrientedData = processSamples(caseData)

dataFrame = pandas.read_csv(xenaFilePath, sep='\t')
myDF = pandas.DataFrame(list(sampleOrientedData.values()))
for col in myDF.columns:
    if col == "samples.submitter_id":
        myDF.rename(columns={col: "sample"}, inplace=True)
        continue
    myDF.rename(columns={col: ".".join((col.split("."))[::-1])}, inplace=True)
myDF["id"] = myDF["case_id"]
if "age_at_diagnosis.diagnoses" in myDF.columns:
    myDF["age_at_earliest_diagnosis.diagnoses.xena_derived"] = myDF['age_at_diagnosis.diagnoses'].apply(lambda x: customMin(x) if isinstance(x, list) else x)
    myDF["age_at_earliest_diagnosis_in_years.diagnoses.xena_derived"] = myDF["age_at_earliest_diagnosis.diagnoses" \
                                                                             ".xena_derived"].apply(lambda x: x/365)

myDF = myDF[list(dataFrame.columns)]

myDF.convert_dtypes().dtypes
dataFrame.convert_dtypes().dtypes

myDF.fillna(pandas.NA, inplace=True)
dataFrame.fillna(pandas.NA, inplace=True)

myDF.sort_values(by=["sample"], inplace=True)
dataFrame.sort_values(by=["sample"], inplace=True)

myDF.reset_index(drop=True, inplace=True)
dataFrame.reset_index(drop=True, inplace=True)

compare()
exit(0)