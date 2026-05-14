#!/usr/bin/env python3
"""Tests for pdif_pwm_scanner.py — PWM construction, scanning, orientation, DC validation, Hall degraded."""

import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import pytest
from PdifFinder.pdif_pwm_scanner import (
    build_pwms,
    scan_sequence,
    validate_dc_sites,
    hall_degraded_scan,
    load_pwms,
    reverse_complement,
)

PDIF1 = "ATTTAACATAAGGGCTGTTATACGAAAT"
XRS38 = "ACTTCGTATAATCGCCATTATGTTAAAT"
SCORE_THRESHOLD = 16.0


@pytest.fixture(scope="module")
def pwms():
    return build_pwms()


@pytest.fixture(scope="module")
def loaded_pwms(pwms):
    return load_pwms()


class TestPWMConstruction:
    """Verify PWM building and JSON loading."""

    def test_build_pwms(self, pwms):
        """build_pwms() returns dict with all required keys, all non-empty."""
        required_keys = ('xerC_pwm', 'xerD_pwm', 'xerC_seeds', 'xerD_seeds')
        for key in required_keys:
            assert key in pwms, f"Missing key: {key}"
            assert pwms[key], f"Empty value for: {key}"

        for pwm_key in ('xerC_pwm', 'xerD_pwm'):
            pwm = pwms[pwm_key]
            for nt in ('A', 'C', 'G', 'T'):
                assert nt in pwm, f"Missing {nt} in {pwm_key}"
                assert len(pwm[nt]) == 11, \
                    f"Expected 11 positions for {pwm_key}[{nt}], got {len(pwm[nt])}"

        assert len(pwms['xerC_seeds']) >= 20
        assert len(pwms['xerD_seeds']) >= 20

    def test_load_pwms_from_json(self, loaded_pwms):
        """load_pwms() returns 5-tuple; seeds >= 46."""
        assert len(loaded_pwms) == 5, \
            f"Expected 5-tuple from load_pwms(), got {len(loaded_pwms)}"

        _, _, _, xerC_seeds, xerD_seeds = loaded_pwms

        assert len(xerC_seeds) >= 46, \
            f"Expected >=46 XerC seeds, got {len(xerC_seeds)}"
        assert len(xerD_seeds) >= 46, \
            f"Expected >=46 XerD seeds, got {len(xerD_seeds)}"


class TestScanSequence:
    """Verify core pdif scanning on known and edge-case sequences."""

    def test_known_xrs38(self, pwms):
        """xrs38 seed: >=1 site, start=0, end=28, raw_orientation='CD'."""
        sites, _ = scan_sequence(
            XRS38,
            pwms['xerC_pwm'], pwms['xerD_pwm'],
            pwms['xerC_seeds'], pwms['xerD_seeds'],
        )
        high = [s for s in sites if s['total_score'] >= SCORE_THRESHOLD]
        assert len(high) >= 1
        s = high[0]
        assert s['start'] == 0
        assert s['end'] == 28
        assert s['raw_orientation'] == 'CD'

    def test_known_pdif1(self, pwms):
        """pdif1: spacer='GGGCTG', CD orientation."""
        sites, _ = scan_sequence(
            PDIF1,
            pwms['xerC_pwm'], pwms['xerD_pwm'],
            pwms['xerC_seeds'], pwms['xerD_seeds'],
        )
        high = [s for s in sites if s['total_score'] >= SCORE_THRESHOLD]
        assert len(high) >= 1
        s = high[0]
        assert s['spacer_seq'] == 'GGGCTG'
        assert s['raw_orientation'] == 'CD'

    def test_flanked_pdif(self, pwms):
        """pdif site flanked by random bases; detected at correct position."""
        prefix = "GATCACAGGTCTAGGATCACAGGTCTAG"
        suffix = "TGCATGCGATCGATGCATGCGATCGA"
        seq = prefix + PDIF1 + suffix

        sites, _ = scan_sequence(
            seq,
            pwms['xerC_pwm'], pwms['xerD_pwm'],
            pwms['xerC_seeds'], pwms['xerD_seeds'],
        )
        high = [s for s in sites if s['total_score'] >= SCORE_THRESHOLD]
        assert len(high) >= 1
        s = high[0]
        assert s['start'] == len(prefix)

    def test_empty_sequence(self, pwms):
        """Empty string yields zero sites."""
        sites, _ = scan_sequence(
            "",
            pwms['xerC_pwm'], pwms['xerD_pwm'],
            pwms['xerC_seeds'], pwms['xerD_seeds'],
        )
        assert len(sites) == 0

    def test_short_sequence(self, pwms):
        """20-mer too short for a full pdif site yields zero sites."""
        sites, _ = scan_sequence(
            "A" * 20,
            pwms['xerC_pwm'], pwms['xerD_pwm'],
            pwms['xerC_seeds'], pwms['xerD_seeds'],
        )
        assert len(sites) == 0

    def test_output_fields(self, pwms):
        """Every output site contains all expected orientation fields."""
        sites, _ = scan_sequence(
            XRS38,
            pwms['xerC_pwm'], pwms['xerD_pwm'],
            pwms['xerC_seeds'], pwms['xerD_seeds'],
        )
        assert len(sites) >= 1

        required_fields = [
            'left_CD_score', 'right_DC_score',
            'orientation_confidence', 'orientation_margin',
            'canonical_orientation', 'raw_orientation',
            'original_orientation', 'motif_supported_orientation',
            'dc_validation_status',
        ]
        for field in required_fields:
            assert field in sites[0], f"Missing required field '{field}'"

    def test_threshold_filter(self, pwms):
        """All returned sites meet the internal 16.5 score threshold."""
        sites, _ = scan_sequence(
            PDIF1,
            pwms['xerC_pwm'], pwms['xerD_pwm'],
            pwms['xerC_seeds'], pwms['xerD_seeds'],
        )
        for s in sites:
            assert s['total_score'] >= 16.5

    def test_reverse_complement_consistency(self, pwms):
        """Reverse complement scan also detects canonical CD sites."""
        sites_fwd, _ = scan_sequence(
            PDIF1,
            pwms['xerC_pwm'], pwms['xerD_pwm'],
            pwms['xerC_seeds'], pwms['xerD_seeds'],
        )
        fwd_cd = [
            s for s in sites_fwd
            if s['canonical_orientation'] == 'CD'
            and s['total_score'] >= SCORE_THRESHOLD
        ]
        assert len(fwd_cd) >= 1

        rc_seq = reverse_complement(PDIF1)
        sites_rev, _ = scan_sequence(
            rc_seq,
            pwms['xerC_pwm'], pwms['xerD_pwm'],
            pwms['xerC_seeds'], pwms['xerD_seeds'],
        )
        rev_cd = [
            s for s in sites_rev
            if s['canonical_orientation'] == 'CD'
            and s['total_score'] >= SCORE_THRESHOLD
        ]
        assert len(rev_cd) >= 1


class TestOrientationScoring:
    """Verify C/D motif orientation scoring fields and values."""

    def test_strong_CD_margin(self, pwms):
        """Known CD site (xrs38) has positive orientation margin and CD label."""
        sites, _ = scan_sequence(
            XRS38,
            pwms['xerC_pwm'], pwms['xerD_pwm'],
            pwms['xerC_seeds'], pwms['xerD_seeds'],
        )
        high = [s for s in sites if s['total_score'] >= SCORE_THRESHOLD]
        assert len(high) >= 1
        s = high[0]
        assert s['orientation_margin'] > 0
        assert s['raw_orientation'] == 'CD'

    def test_cd_dc_scores_present(self, pwms):
        """left_CD_score and right_DC_score are integers."""
        sites, _ = scan_sequence(
            XRS38,
            pwms['xerC_pwm'], pwms['xerD_pwm'],
            pwms['xerC_seeds'], pwms['xerD_seeds'],
        )
        assert len(sites) >= 1
        s = sites[0]
        assert isinstance(s['left_CD_score'], int)
        assert isinstance(s['right_DC_score'], int)


class TestDCValidation:
    """Verify validate_dc_sites() status strings under various scenarios."""

    def test_validate_dc_supported(self):
        """motif=DC + GB=DC -> validated_DC."""
        site = {
            'start': 0, 'total_score': 20.0,
            'motif_supported_orientation': 'DC', 'orientation_margin': 5.0,
            'orientation_reason': '', 'dc_validation_status': 'not_GB_DC',
        }
        result = validate_dc_sites([site], {0: 'DC'})
        assert result[0]['dc_validation_status'] == 'validated_DC'

    def test_validate_dc_rejected(self):
        """motif=CD + GB=DC -> GB_DC_not_supported_by_motif."""
        site = {
            'start': 0, 'total_score': 20.0,
            'motif_supported_orientation': 'CD', 'orientation_margin': 3.0,
            'orientation_reason': '', 'dc_validation_status': 'not_GB_DC',
        }
        result = validate_dc_sites([site], {0: 'DC'})
        assert result[0]['dc_validation_status'] == 'GB_DC_not_supported_by_motif'

    def test_validate_ambiguous(self):
        """motif=ambiguous + GB=DC -> ambiguous_DC."""
        site = {
            'start': 0, 'total_score': 20.0,
            'motif_supported_orientation': 'ambiguous', 'orientation_margin': 0,
            'orientation_reason': 'equal_scores', 'dc_validation_status': 'not_GB_DC',
        }
        result = validate_dc_sites([site], {0: 'DC'})
        assert result[0]['dc_validation_status'] == 'ambiguous_DC'

    def test_non_dc_unchanged(self):
        """GB=CD sites remain not_GB_DC."""
        site = {
            'start': 0, 'total_score': 20.0,
            'motif_supported_orientation': 'CD', 'orientation_margin': 5.0,
            'orientation_reason': '', 'dc_validation_status': 'not_GB_DC',
        }
        result = validate_dc_sites([site], {0: 'CD'})
        assert result[0]['dc_validation_status'] == 'not_GB_DC'


class TestHallDegraded:
    """Verify hall_degraded_scan() guided and unguided modes."""

    def test_guided_search_detects_known(self, pwms):
        """Known pdif site found in 'high' or 'medium' with expected positions."""
        result = hall_degraded_scan(
            PDIF1,
            pwms['xerC_pwm'], pwms['xerD_pwm'],
            pwms['xerC_seeds'], pwms['xerD_seeds'],
            expected_positions=[0],
        )
        assert len(result['high']) >= 1 or len(result['medium']) >= 1

    def test_unguided_scan(self, pwms):
        """Without expected positions, hall_degraded_scan still finds sites."""
        result = hall_degraded_scan(
            PDIF1,
            pwms['xerC_pwm'], pwms['xerD_pwm'],
            pwms['xerC_seeds'], pwms['xerD_seeds'],
        )
        assert len(result['high']) >= 1
