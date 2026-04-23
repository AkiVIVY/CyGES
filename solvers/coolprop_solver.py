import CoolProp as CP

from core.units import jpkg_to_kjpkg, kjpkg_to_jpkg, kpa_to_pa, pa_to_kpa


class CoolPropSolver:
    """
    新架构原生物性求解器（替代对旧propertySolver.py的直接依赖）。
    对外统一使用工程单位:
    - T: K
    - P: kPa
    - H: kJ/kg
    - S: kJ/(kg.K)
    """

    _PAIR_TO_INPUT = {
        "HP": CP.HmassP_INPUTS,
        "TP": CP.PT_INPUTS,
        "HS": CP.HmassSmass_INPUTS,
        "PS": CP.PSmass_INPUTS,
    }

    def __init__(self, substance_list: list[str]):
        self.substance_list = substance_list
        self._states = {}
        for substance in substance_list:
            try:
                self._states[substance] = CP.AbstractState("HEOS", substance)
            except Exception as exc:
                raise ValueError(f"invalid CoolProp substance: {substance}") from exc

    def solve(self, pair: str, substance: str, input1: float, input2: float) -> dict:
        pair = pair.upper()
        if pair not in self._PAIR_TO_INPUT:
            raise ValueError(f"unsupported input pair: {pair}")
        if substance not in self._states:
            raise ValueError(f"substance not initialized: {substance}")

        state = self._states[substance]
        cp_pair = self._PAIR_TO_INPUT[pair]
        v1, v2 = self._to_coolprop_inputs(pair, input1, input2)
        state.update(cp_pair, v1, v2)

        return {
            "T": state.T(),
            "P": pa_to_kpa(state.p()),
            "H": jpkg_to_kjpkg(state.hmass()),
            "S": jpkg_to_kjpkg(state.smass()),
            "substance": substance,
        }

    @staticmethod
    def _to_coolprop_inputs(pair: str, input1: float, input2: float) -> tuple[float, float]:
        if pair == "HP":
            return kjpkg_to_jpkg(input1), kpa_to_pa(input2)
        if pair == "TP":
            return kpa_to_pa(input2), input1
        if pair == "HS":
            return kjpkg_to_jpkg(input1), kjpkg_to_jpkg(input2)
        if pair == "PS":
            return kpa_to_pa(input1), kjpkg_to_jpkg(input2)
        raise ValueError(f"unsupported input pair: {pair}")

