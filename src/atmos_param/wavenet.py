import os
import pickle
import time
import traceback
from typing import Any, Dict, Union

import numpy as np
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers

"""
Choose input features from union of:
------------------------------------
# WaveNet ["temp", "hght", "ucomp", "vcomp", "omega", "slp", "lat", "lon"]
# MiMA ["lat", "lon", "lat", "pfull", "zfull", "temp", "uuu", "vvv", "time", "delt"]
"""
features_gwfu = ["uuu", "temp"]
features_gwfv = ["vvv", "temp"]
features_combined = ["temp", "uuu", "vvv", "lat", "lon"]
model_path = '/home/zespinos/models/wavenet_2/' # TODO: use '/home/zespinos/models/wind_temp'
ZERO_LAYERS = 7


def build_input_tensors(model: str, args: Dict[str, np.ndarray]):
    if model == 'gwfu':
        features = features_gwfu
        output_scaler = scaler_gwfu
    elif model == "combined":
        features = features_combined
        output_scaler = scaler_combined
    else:
        features = features_gwfv
        output_scaler = scaler_gwfv

    tensors = []
    for feat in features:
        if isinstance(args[feat], int):
            feat = np.full((args["temp"].shape[0]*args["temp"].shape[1], 1), fill_value=args[feat])
        else:
            feat = args[feat].reshape(args[feat].shape[0]*args[feat].shape[1], args[feat].shape[2])
        tensors.append(feat)

    tensors = np.concatenate(tensors, axis=1)
    tensors = scaler_tensors.transform(tensors)

    return tensors


def from_pickle(path: Union[str, os.PathLike]) -> Any:
    try:
        with open(path, 'rb') as handle:
            b = pickle.load(handle)
    except:
        traceback.print_exc()
    else:
        return b


def generate_predictions(wavenet, scaler, tensors, args):
    # TODO: Adding zeros to end (i.e. assuming 0 index is top [.18hPa] and 33 is bottom [475hPa])
    predictions = wavenet.predict(tensors)
    predictions = scaler.inverse_transform(predictions)

    predictions = np.concatenate(
        (
            predictions[:, :33],
            np.zeros((predictions.shape[0], ZERO_LAYERS)),
            predictions[:,33:],
            np.zeros((predictions.shape[0], ZERO_LAYERS))
        ), axis=1
    )

    predictions = predictions.reshape((args[6].shape[0], args[6].shape[1], args[6].shape[2]*2))
    predictions = np.asfortranarray(predictions)

    return (predictions[:,:,:40], predictions[:,:,40:])


def predict(*args):
    # Build Input Vectors
    tensors = build_input_tensors(
        model="combined",
        args={
            'temp': args[5],  # temp
            'zfull': args[4], # zfull
            'uuu': args[6],   # uuu
            'vvv': args[7],   # vvv
            'lat': args[0],   # lat
            'lon': args[1],   # lon
        }
    )
    predictions = generate_predictions(
        wavenet=wavenet_combined,
        scaler=scaler_combined,
        tensors=tensors,
        args=args,
    )

    return predictions


def init_wavenet():
    # Load Model Once
    global wavenet_combined
    global scaler_combined
    global scaler_tensors

    wavenet_combined = tf.keras.models.load_model(os.path.join(model_path, 'wavenet_2.hdf5'), compile=False)
    wavenet_combined.compile(loss='logcosh', optimizer="adam")
    # Load Scalars Once
    scaler_combined = from_pickle(os.path.join(model_path, 'combined_scaler.pkl'))
    scaler_tensors = from_pickle(os.path.join(model_path, 'tensors_scaler.pkl'))


###### Development Testing Only #####
def local_testing():
    dd_args = from_pickle(os.path.join('./dd_args', 'args.pkl'))
    gwfu, gwfv = predict(*list(dd_args.values()))
    print("GWFU shape: ", gwfu.shape)
    print("GWFV shape: ", gwfv.shape)



#####################################
#init_wavenet()
#local_testing()
