#!/usr/bin/env python3 -u
# Copyright (c) Facebook, Inc. and its affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

import argparse
import pathlib

import torch

from esm import Alphabet, FastaBatchedDataset, ProteinBertModel, pretrained



def main(model_location,input_file,input_embedding,repr_layer):
    #repr_layer=[*range(0,34,1)]
    nogpu='--nogpu'
    toks_per_batch=4096
    include=['mean','per_tok']
    model, alphabet = pretrained.load_model_and_alphabet(model_location)
    model.eval()
    if torch.cuda.is_available() and not nogpu:
        model = model.cuda()
        print("Transferred model to GPU")

    dataset = FastaBatchedDataset.from_file(input_file)
    batches = dataset.get_batch_indices(int(toks_per_batch), extra_toks_per_seq=1)
    data_loader = torch.utils.data.DataLoader(
        dataset, collate_fn=alphabet.get_batch_converter(), batch_sampler=batches
    )
    print(f"Read {input_file} with {len(dataset)} sequences")

    return_contacts = "contacts" in include

    assert all(-(model.num_layers + 1) <= i <= model.num_layers for i in repr_layer)
    repr_layers = [(i + model.num_layers + 1) % (model.num_layers + 1) for i in repr_layer]

    with torch.no_grad():
        for batch_idx, (labels, strs, toks) in enumerate(data_loader):
            print(
                f"Processing {batch_idx + 1} of {len(batches)} batches ({toks.size(0)} sequences)"
            )
            if torch.cuda.is_available() and not nogpu:
                toks = toks.to(device="cuda", non_blocking=True)

            # The model is trained on truncated sequences and passing longer ones in at
            # infernce will cause an error. See https://github.com/facebookresearch/esm/issues/21
            
            out = model(toks, repr_layers=repr_layers, return_contacts=return_contacts)

            logits = out["logits"].to(device="cpu")
            representations = {
                layer: t.to(device="cpu") for layer, t in out["representations"].items()
            }
            if return_contacts:
                contacts = out["contacts"].to(device="cpu")

            for i, label in enumerate(labels):
                output_file = input_embedding + f"{label}.pt"
                result = {"label": label}
                # Call clone on tensors to ensure tensors are not views into a larger representation
                # See https://github.com/pytorch/pytorch/issues/1995
                if "per_tok" in include:
                    result["representations"] = {
                        layer: t[i, 1 : len(strs[i]) + 1].clone()
                        for layer, t in representations.items()
                    }
                if "mean" in include:
                    result["mean_representations"] = {
                        layer: t[i, 1 : len(strs[i]) + 1].mean(0).clone()
                        for layer, t in representations.items()
                    }
                if "bos" in include:
                    result["bos_representations"] = {
                        layer: t[i, 0].clone() for layer, t in representations.items()
                    }
                if return_contacts:
                    result["contacts"] = contacts[i, : len(strs[i]), : len(strs[i])].clone()

                torch.save(
                    result,
                    output_file,
                )


