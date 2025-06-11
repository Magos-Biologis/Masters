#!../.venv/bin/python3
import json
import os
import re

import numpy as np
import pandas as pd
from pylab import *

# figure_env = str(os.getenv("THESIS_FIGURE_PATH"))
fig_env = str(os.getenv("THESIS_FIGURE_PATH"))
fpe_env = str(os.getenv("FPE_FIGURE_ENV"))
phs_env = str(os.getenv("PHS_FIGURE_ENV"))
ode_env = str(os.getenv("ODE_FIGURE_ENV"))

data_env = str(os.getenv("THESIS_DATA_PATH"))

ode_dir = os.path.join(fpe_env, "ode")
fpe_dir = os.path.join(fpe_env, "five_var")


## Taken from my parameter flag in the sim script
filter_captures = re.compile(r"(\w[^=]*)=\s?([^ ,]*)")
# kwarg_string_parser = re.compile(r"(\w[^=]*)=\s?([^ ,]*)")


def parse_filters(string):
    matches: list[tuple] = filter_captures.findall(string)
    parameters = [
        (str(key).replace("-", "_").replace("'", "p").replace(" ", ""), float(value))
        for key, value in matches
    ]
    return parameters


# @njit
def string_or_float(input: str) -> str | float:
    if input.lower() == "true" or input.lower() == "false":
        return bool(input)
    else:
        try:
            return float(input)
        except:
            return str(input)


def parse_kwarg_string(string):
    matches: list[tuple] = filter_captures.findall(string)
    parameters = [
        (str(key).replace(" ", ""), string_or_float(value)) for key, value in matches
    ]
    return dict(parameters)

    choices = (
        [
            "markov",
            "ode",
            "phase",
            "ssa",
        ],
    )


## Compiling the defaults and the choice of parameters


proper_name_capture: re.Pattern = re.compile(r"((ODE)|(Novel)|(Simple)|(Three)|(Ill))")
metadata_captures: re.Pattern = re.compile(r"(\w[^=]*)=([^_]*)_")
ratio_pattern: re.Pattern = re.compile(r"R(?P<ratio>b.n)R")
capture_patterns: re.Pattern = re.compile(
    r"(?P<data_source>^[^M]*)"  # phase, ssa, etc.
    r"M(?P<model_name>[^P]*)"  # Model specifier
    r"P(?P<rawmetadata>[^SR]*)"  # Parameters and
    r"[^S]*"  # a stupid way to ignore the ratio
    r"S(?P<steps_taken>[^T]*)"  # Step count
    r"I(?P<initial_condition>[^CT]*)C?"  # Initial conditions
    r"T(?P<date>\d+)"
    r"\.(?P<filetype>\w*)$"
)


def parse_metadata(string):
    matches: list[tuple] = metadata_captures.findall(string)
    parameters = [
        (str(key).replace("-", "_").replace("km", "k_"), float(value))
        for key, value in matches
    ]
    return dict(parameters)


data_source_map = {
    "ssa": "stochastic simulation algorithm",
    "phase": "model phase space",
    "ode": "model integration",
    "markov": "markov chain simulation",
}

model_map = {
    "jesper": "ThreeCompartment",
    "2L": "SimpleChemical1Par2Var",
    "2S": "BrokenSimpleChemical2Par2Var",
    "ode_2_2": "ODE2Par2Var",
    "ode_3_2": "ODE3Par2Var",
    "ode_3_3": "ODE3Par3Var",
    "ode_5_2": "ODE5Par2Var",
    "ode_5_3": "ODE5Par3Var",
    "ode_7_3": "ODE7Par2Var",
    "ode_8_3": "ODE8Par3Var",
    "5_2": "Novel5Par2Var",
    "5_3": "Novel5Par3Var",
}


records = []
file_names = os.listdir(data_env)
for fn in file_names:
    correct_format = proper_name_capture.search(fn)
    if correct_format is not None:
        # print("Already Updated Filename")
        continue
        old_name = os.path.join(data_env, fn)

        ## This monstrocity is an attempt at loading correctly encoded json
        ## data into a python `dict()`, its important because otherwise
        ## I have no way to access the encoded data, as some of it is
        ## corrupted
        try:
            numpy_obj = np.load(old_name)
        except:
            pass
        else:
            try:
                numpy_json_obj = numpy_obj["metadata"]
            except KeyError:
                continue
            else:
                if type(numpy_json_obj) is np.ndarray:
                    numpy_bytes = numpy_json_obj.tobytes(order="C")
                    try:
                        json_obj = json.loads(numpy_bytes)
                    except:
                        print("json didn't approve")
                    else:
                        pass
                        # print("json likey")
                        # print(numpy_bytes)
                        # print(len(str(json_obj["date"])))
        continue

    rec = dict()

    matches = capture_patterns.match(fn)
    ratio = ratio_pattern.search(fn)

    if matches is not None:
        rec.update(matches.groupdict())

    if ratio is not None:
        rec.update(ratio.groupdict())

    rec["original_file_name"] = fn
    new_file_name = rec["data_source"]

    try:
        rec["data_source"] = data_source_map[rec["data_source"]]
    except KeyError:
        rec.update(data_source=None)

    try:
        rec["model_name"] = model_map[rec["model_name"]]
    except KeyError:
        rec["model_name"] = "IllDefinedModel"

    new_file_name += "M{}".format(rec["model_name"])
    new_file_name += "L{}".format("python")

    ## Making the metadata 'real' entries
    ## maybe I should make it into a sub-dict?
    ## Nah, I can always do that later
    try:
        parameter_data = parse_metadata(rec["rawmetadata"])
    except KeyError:
        rec.update(parameters=None)
    else:
        rec.update(parameters=parameter_data)

    ## Ensuring the right datatype for the step count
    try:
        rec["steps_taken"] = int(float(rec["steps_taken"].replace("S", "")))
    except KeyError:
        rec.update(steps=None)
    except ValueError:
        rec.update(steps=None)

    ## Initial condition and number of variable handling
    ## The idea is that if I have a well defined initial condition, I also
    ## get the number of variables from the length
    try:
        rec["initial_condition"] = [
            int(v)
            for v in rec["initial_condition"].replace("[", "").replace("]", "").split()
        ]
    except KeyError:
        rec.update(initial_condition=None)
    else:
        rec.update(number_of_variables=int(len(rec["initial_condition"])))

    new_file_name += "T{}".format(rec["date"])
    ## Consistency of time formating
    if len(rec["date"]) > 10:
        rec["date"] = float(rec["date"][:10] + "." + rec["date"][10:])
    else:
        try:
            rec["date"] = int(rec["date"])
        except:
            print("poop")
            exit()

    ## Number of particle handling
    try:
        quantity = int(rec["parameters"].pop("num"))
    except KeyError:
        rec.update(number_of_particles=int(sum(rec["initial_condition"])))
    except ValueError:
        rec.update(number_of_particles=None)
    else:
        rec.update(number_of_particles=quantity)

    ## This is just dumping some needless details
    ## The pop function removing the value if it exists, returning a default if not
    _ = rec.pop("filetype", "")
    _ = rec["parameters"].pop("k_21.0_k3", "")
    _ = rec["parameters"].pop("k_31.0_k4", "")
    _ = rec["parameters"].pop("k_41.0_k5", "")

    rec.update(run=1)

    ## Actually defining the json file, after being confident I have all the
    ## needed metadata, or in the very least that I can construct anything
    ## missing.
    ## Ensuring its encoded in utf-8 so as to maintain consistency with
    ## the Julia output

    metadata_json = json.dumps(rec, ensure_ascii=False).encode("utf-8")
    new_file_name += ".npz"

    old_name = os.path.join(data_env, fn)
    new_name = os.path.join(data_env, new_file_name)

    try:
        original_data = np.load(old_name)
    except:
        print("no data")
        continue

    try:
        np.savez(
            old_name,
            time=original_data["time"],
            states=original_data["states"],
            metadata=metadata_json,
        )
        # new_data = np.load(old_name)
        # json_data = json.loads(new_data["metadata"].tobytes(order="C"))
    except:
        try:
            # print(original_data)
            # print(original_data["c1"])
            # continue
            np.savez(
                old_name,
                time=original_data["time"],
                solution=original_data["solutions"],
                metadata=metadata_json,
            )
        except:
            ## I forgot that I set the time to 0 for all the phase stuff,
            ## lowkey just erased all that data
            ## Whatever, it works fine cause I don't need to worry about
            ## any randomness
            try:
                for_npz = {
                    "c1": original_data["c1"],
                    "c2": original_data["c2"],
                    "dU": original_data["dU"],
                    "dV": original_data["dV"],
                }
            except:
                print("Failed to load data")
                continue
            else:
                try:
                    nulls = {
                        "c1_nullcline": original_data["c1_nullcline"],
                        "c2_nullcline": original_data["c2_nullcline"],
                    }
                except:
                    pass
                else:
                    for_npz.update(nulls)

                np.savez(
                    old_name,
                    **for_npz,
                    metadata=metadata_json,
                )

    # print(old_name)
    # print(new_name)

    try:
        os.rename(old_name, new_name)
    except:
        print("Failed to rename file")
        continue

print("Done Renaming")
exit()


# no_ratio = raw_frame["ratio"].isna()
# raw_frame.loc[no_ratio, "ratio"] = ""


### Now we can mutate the dictionaries in memory and dynamically create columns
defined_metadata = raw_frame["metadata"].map(lambda d: len(d) > 0)
raw_frame = raw_frame.loc[defined_metadata]

## That dict.pop has a default for missing entries is a literal godsend
ssa_index = raw_frame.loc[raw_frame["datasource"] == "ssa"].index
raw_frame["count"] = [d.pop("num", np.nan) for d in raw_frame.loc[:, "metadata"]]


## From the list, we then count the number of times each model shows up,
## assigning the number of each step accordingly.
## From which multi-index is generated for the dataframe
## so as to have an easier time choosing which file to plot
raw_frame["model_index"] = raw_frame.groupby("model").cumcount()
file_frame = raw_frame.set_index(["datasource", "model", "model_index"]).sort_index()

parameters = file_frame["metadata"].apply(pd.Series)
file_frame = pd.concat([file_frame, parameters], axis=1)
sourced_frame = (
    file_frame.loc[(data_source, model)]
    .dropna(axis="columns", how="all")
    .reset_index(drop=True)
)


def renaming_function(entry):
    pass


#######################
#   File name stuff   #
#######################


## The base name
# file_name: str = "{}:{}:".format(data_source, model).replace("_", "-")
# if model == "5_2":
#     para_version = data_frame["ratio"]
#     file_name += "{}:".format(para_version).replace("$", "")
#
# # Format rapidfire
# file_name += "I{}:".format(m_index)
# if args.plot_fixedpoints:
#     file_name += "C{}".format(c_data)
#     file_name += "{}:".format("fp")
# else:
#     file_name += "C{}:".format(c_data)
#
# file_name = data_source

##########################


## So we assign it to a dataframe for easy access
# raw_frame = pd.DataFrame(records)
# raw_frame.sort_values(by="t", inplace=True)  ## Sorting by unix epoch time, ascending


##########################

exit()

print("Done sorting")
exit()
