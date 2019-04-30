from .code_info import features


def has_features(*args):
    """Tests whether a list of features is a subset of the compiled-in features"""

    if len(args) == 1 and not isinstance(args[0], str) and hasattr(args[0], "__iter__"):
        return set(args[0]) < set(features())

    return set(args) < set(features())


def missing_features(*args):
    """Returns a list of the missing features in the argument"""

    if len(args) == 1 and not isinstance(args[0], str) and hasattr(args[0], "__iter__"):
            return set(args[0]) - set(features())

    return set(args) - set(features())


def assert_features(*args):
    """Raises an exception when a list of features is not a subset of the compiled-in features"""

    if not has_features(*args):
        raise Exception(
            "Missing features " + ", ".join(missing_features(*args)))
