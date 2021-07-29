import pint
from equilibrator_api import (
    ComponentContribution,
    Reaction,
    default_physiological_ionic_strength,
    default_physiological_p_h,
    default_physiological_p_mg,
    default_physiological_temperature,
)
from equilibrator_api.phased_reaction import PhasedReaction
from equilibrator_assets.compounds import Compound
from equilibrator_assets.local_compound_cache import LocalCompoundCache
from equilibrator_cache.compound_cache import CompoundCache
from pymongo import MongoClient
from sqlalchemy import create_engine
from typing import Union


class MINE_thermo:
    def __init__(
        self,
        mongo_uri: Union[str, None] = None,
        postgres_uri: str = "postgresql:///eq_compounds",
    ):
        self.CC = ComponentContribution()

        if mongo_uri:
            self.mongo_uri = mongo_uri
            self.client = MongoClient(mongo_uri)
        else:
            self.mongo_uri = "localhost:27017"
            self.client = MongoClient()
        self.core = self.client.core

        self.postres_uri = postgres_uri
        self.pc = LocalCompoundCache()
        self.pc.ccache = CompoundCache(create_engine(postgres_uri))

        self.water = self.pc.get_compounds("O")

    def get_eQ_compound_from_cid(self, c_id: str) -> Compound:
        """Get an equilibrator compound for a given c_id.

        Parameters
        ----------
        c_id : str
            compound ID for MongoDB lookup of a compound.

        Returns
        -------
        equilibrator_assets.compounds.Compound
            eQuilibrator Compound
        """
        compound = self.core.compounds.find_one({"_id": c_id}, {"SMILES": 1})
        compound_smiles = compound["SMILES"]
        eQ_compound = self.pc.get_compounds(
            compound_smiles, bypass_chemaxon=True, save_empty_compounds=True
        )


        return eQ_compound

    def standard_dg_formation_from_cid(self, c_id: str) -> pint.Measurement:
        """Get standard ∆Gf for a compound.

        Parameters
        ----------
        c_id : str
            Compound ID to get the ∆Gf for.

        Returns
        -------
        pint.Measurement
            ∆Gf for a compound.
        """
        eQ_cpd = self.get_eQ_compound_from_cid(c_id)
        dgf = self.CC.standard_dg_formation(eQ_cpd)

        return dgf

    def get_eQ_reaction_from_rid(self, r_id: str, mine: str) -> PhasedReaction:
        """Get an eQuilibrator reaction object from an r_id.

        Parameters
        ----------
        r_id : str
            Reaction id to get object for.
        mine : str
            Database to look for reaction in.

        Returns
        -------
        PhasedReaction
            eQuilibrator reactiono to calculate ∆Gr with.
        """
        mine = self.client[mine]
        reaction_info = mine.reactions.find_one({"_id": r_id})
        reactants = reaction_info["Reactants"]
        products = reaction_info["Products"]

        lhs = " + ".join(f"{r[0]} {r[1]}" for r in reactants)
        rhs = " + ".join(f"{p[0]} {p[1]}" for p in products)
        reaction_string = " => ".join([lhs, rhs])

        compounds = set([r[1] for r in reactants])
        compounds.update(tuple(p[1] for p in products))

        eQ_compound_dict = {
            c_id: self.get_eQ_compound_from_cid(c_id) for c_id in compounds
        }

        if "X73bc8ef21db580aefe4dbc0af17d4013961d9d17" not in compounds:
            eQ_compound_dict["water"] = self.water

        eq_reaction = Reaction.parse_formula(eQ_compound_dict.get, reaction_string)

        return eq_reaction

    def physiological_dg_prime_from_rid(self, r_id: str, mine: str):
        """Calculate the ∆G'physiological of a reaction.

        Parameters
        ----------
        r_id : str
            ID of the reaction to calculate.
        mine : str
            MINE the reaction is found in.

        Returns
        -------
        pint.Measurement
            The calculated ∆G'physiological.
        """
        eQ_reactiton = self.get_eQ_reaction_from_rid(r_id, mine)
        dgrp = self.CC.physiological_dg_prime(eQ_reactiton)

        return dgrp

    def standard_dg_from_rid(self, r_id: str, mine: str):
        """Calculate the ∆Go of a reaction.

        Parameters
        ----------
        r_id : str
            ID of the reaction to calculate.
        mine : str
            MINE the reaction is found in.

        Returns
        -------
        pint.Measurement
            The calculated ∆Go.
        """
        eQ_reactiton = self.get_eQ_reaction_from_rid(r_id, mine)
        dgrp = self.CC.standard_dg_prime(eQ_reactiton)

        return dgrp

    # Keep this just in case we decided to modify the CC at some point.
    # def _reset_CC_conditions(self) -> None:
    #     self.CC.p_h = default_physiological_p_h
    #     self.CC.p_mg = default_physiological_p_mg
    #     self.CC.temperature = default_physiological_temperature
    #     self.CC.ionic_strength = default_physiological_ionic_strength