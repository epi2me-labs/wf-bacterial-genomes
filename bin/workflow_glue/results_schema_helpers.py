"""Helpers for the auto-generated schema code."""
from pydantic import BaseModel as PydanticBaseModel


class BaseModel(PydanticBaseModel):
    """Extend base model."""

    class Config:
        """Config items for the pydantic code."""

        # make enums json serializable
        use_enum_values = True
