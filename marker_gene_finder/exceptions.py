"""Documented exceptions raised by the public v2 API."""


class MarkerFinderError(Exception):
    """Base class for public Marker Gene Finder failures."""


class ConfigurationError(MarkerFinderError, ValueError):
    """An analysis configuration is invalid."""


class InputValidationError(MarkerFinderError, ValueError):
    """An input dataset cannot be analyzed with the requested configuration."""


class SchemaError(MarkerFinderError, ValueError):
    """A serialized v2 document has an unsupported or invalid schema."""
