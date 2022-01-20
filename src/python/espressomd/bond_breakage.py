from .script_interface import script_interface_register, ScriptObjectMap, ScriptInterfaceHelper


@script_interface_register
class BreakageSpec(ScriptInterfaceHelper):
    _so_name = "BondBreakage::BreakageSpec"


class BreakageSpecs(ScriptObjectMap):
    _so_name = "BondBreakage::BreakageSpecs"
    _key_type = int
