def validate_system_spec(spec: dict) -> None:
    required_top_level = ["hot_streams", "cold_streams", "closed_cycles"]
    for key in required_top_level:
        if key not in spec:
            raise ValueError(f"missing required key: {key}")

    for item in spec["hot_streams"] + spec["cold_streams"]:
        for key in ["stream_id", "fluid", "m_dot"]:
            if key not in item:
                raise ValueError(f"stream missing required key: {key}")

    for item in spec["closed_cycles"]:
        for key in ["cycle_id", "fluid", "boundary", "levels"]:
            if key not in item:
                raise ValueError(f"closed_cycle missing required key: {key}")
        levels = item["levels"]
        for level_key in ["TLevel", "PLevel"]:
            if level_key not in levels:
                raise ValueError(f"closed_cycle levels missing required key: {level_key}")
            if not isinstance(levels[level_key], (list, tuple)) or len(levels[level_key]) < 1:
                raise ValueError(f"{level_key} must contain at least one value")

