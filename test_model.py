import glob

from model import Syntax, Part
from pydantic import ValidationError
import unittest


class TestSyntax(unittest.TestCase):
    def test_existing_syntax(self):
        for syntax_file in glob.glob("syntaxes/*/syntax.json"):
            with open(syntax_file, "r") as f:
                Syntax.model_validate_json(f.read())

    def _get_valid_syntax_dict(self):
        """Helper method to create a minimal valid syntax dictionary."""
        return {
            "syntaxName": "Test Syntax",
            "assemblyEnzyme": "BsaI",
            "domesticationEnzyme": "BsmBI",
            "relatedDois": ["10.1000/xyz123"],
            "submitters": ["0000-0000-0000-0000"],
            "overhangNames": {"ACGT": "test_overhang"},
            "parts": [
                {
                    "id": 1,
                    "name": "part1",
                    "left_overhang": "ACGT",
                    "right_overhang": "CGTA",
                }
            ],
        }

    def test_validate_submitters(self):
        syntax_dict = self._get_valid_syntax_dict()
        Syntax.model_validate(syntax_dict)

        # Invalid ORCID formats
        self.assertRaises(
            ValidationError,
            Syntax.model_validate,
            {**syntax_dict, "submitters": ["0000-0000-0000"]},
        )
        self.assertRaises(
            ValidationError,
            Syntax.model_validate,
            {**syntax_dict, "submitters": ["invalid"]},
        )

    def test_validate_overhang_names(self):
        syntax_dict = self._get_valid_syntax_dict()
        Syntax.model_validate(syntax_dict)

        # Invalid overhang names (not 4 characters)
        self.assertRaises(
            ValidationError,
            Syntax.model_validate,
            {**syntax_dict, "overhangNames": {"ACG": "test"}},
        )
        self.assertRaises(
            ValidationError,
            Syntax.model_validate,
            {**syntax_dict, "overhangNames": {"ACGTG": "test"}},
        )

    def test_validate_parts_unique_ids(self):
        syntax_dict = self._get_valid_syntax_dict()
        Syntax.model_validate(syntax_dict)

        # Invalid parts with duplicate IDs
        self.assertRaises(
            ValidationError,
            Syntax.model_validate,
            {
                **syntax_dict,
                "parts": [
                    {
                        "id": 1,
                        "name": "part1",
                        "left_overhang": "ACGT",
                        "right_overhang": "CGTA",
                    },
                    {
                        "id": 1,
                        "name": "part2",
                        "left_overhang": "TTTT",
                        "right_overhang": "AAAA",
                    },
                ],
            },
        )


class TestPart(unittest.TestCase):
    def test_part_validation(self):
        part_dict = {
            "id": 1,
            "name": "1",
            "info": "Assembly connector",
            "glyph": "three-prime-sticky-restriction-site",
            "left_overhang": "CCCT",
            "right_overhang": "AACG",
            "left_inside": "A",
            "right_inside": "ATTTTTTTT",
            "left_codon_start": 0,
            "right_codon_start": 0,
            "color": "#84c5de"
        }
        Part.model_validate(part_dict)

        # Invalid overhangs or inside
        for side in ["left", "right"]:
            self.assertRaises(ValidationError, Part.model_validate, {**part_dict, f"{side}_overhang": "CCCTT"})
            self.assertRaises(ValidationError, Part.model_validate, {**part_dict, f"{side}_overhang": "NNNN"})
            self.assertRaises(ValidationError, Part.model_validate, {**part_dict, f"{side}_overhang": ""})
            self.assertRaises(ValidationError, Part.model_validate, {**part_dict, f"{side}_inside": "NNNN"})

            Part.model_validate({**part_dict, f"{side}_inside": "ATTTTTTTT"})
            Part.model_validate({**part_dict, f"{side}_inside": ""})
            Part.model_validate({**part_dict, f"{side}_codon_start": 0})
            Part.model_validate({**part_dict, f"{side}_codon_start": 1000})

            # Invalid codon start
            self.assertRaises(ValidationError, Part.model_validate, {**part_dict, f"{side}_codon_start": -1})

        # Invalid color
        self.assertRaises(ValidationError, Part.model_validate, {**part_dict, "color": "invalid"})

    def test_default_values(self):
        part = Part(
            id=1,
            name="1",
            left_overhang="CCCT",
            right_overhang="AACG",
        )
        self.assertEqual(part.info, "")
        self.assertEqual(part.glyph, "")
        self.assertEqual(part.left_inside, "")
        self.assertEqual(part.right_inside, "")
        self.assertEqual(part.left_codon_start, 0)
        self.assertEqual(part.right_codon_start, 0)
        self.assertEqual(part.color, "")
