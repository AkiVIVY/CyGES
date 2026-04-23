"""
项目统一工程单位:
- T: K
- P: kPa
- H: kJ/kg
- S: kJ/(kg.K)
"""


def pa_to_kpa(value_pa: float) -> float:
    return value_pa / 1000.0


def kpa_to_pa(value_kpa: float) -> float:
    return value_kpa * 1000.0


def jpkg_to_kjpkg(value_jpkg: float) -> float:
    return value_jpkg / 1000.0


def kjpkg_to_jpkg(value_kjpkg: float) -> float:
    return value_kjpkg * 1000.0

