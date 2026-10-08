"""
tests for pw2cw/postcw command line interface and error messages.
"""

import pytest
from click.testing import CliRunner

from symclosestwannier.__init__ import __version__
from symclosestwannier.scripts.pw2cw import cmd as pw2cw
from symclosestwannier.scripts.postcw import cmd as postcw
from symclosestwannier.util.exceptions import SymCWFileNotFoundError, SymCWInputError

COMMANDS = [pw2cw, postcw]


# ==================================================
@pytest.mark.parametrize("command", COMMANDS)
def test_version(command, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    result = CliRunner().invoke(command, ["-v", "seed"])

    assert result.exit_code == 0
    assert result.output.strip() == f"SymClosestWannier: {__version__}"


# ==================================================
@pytest.mark.parametrize("command", COMMANDS)
def test_input_format(command):
    result = CliRunner().invoke(command, ["-i"])

    assert result.exit_code == 0
    assert "outdir" in result.output


# ==================================================
@pytest.mark.parametrize("command", COMMANDS)
def test_missing_seedname(command):
    result = CliRunner().invoke(command, [])

    assert result.exit_code == 2
    assert "missing SEEDNAME" in result.output


# ==================================================
@pytest.mark.parametrize("command", COMMANDS)
def test_missing_cwin(command, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    result = CliRunner().invoke(command, ["seed.cwin"])

    assert result.exit_code == 1
    assert "Error: cannot find the cwin file:" in result.output
    assert "seed.cwin" in result.output
    assert "Traceback" not in result.output


# ==================================================
@pytest.mark.parametrize(
    "extra, message",
    [
        ("delta = 1e-12", "invalid keyword = delta was given. did you mean 'cwf_delta'"),
        ("ket_amn = [pz@C_1,pz@C_2]", "invalid ket_amn = [pz@C_1,pz@C_2] was given. each ket must be"),
        ("proj_min = abc", "invalid value proj_min = abc"),
        ("restart = foo", "invalid restart = foo was given."),
        ("fermi_energy = 0.0", "invalid keyword = fermi_energy was given. fermi_energy must be given in seedname.win."),
    ],
)
def test_invalid_cwin(make_case, extra, message):
    workdir = make_case("ch4_sl", extra=extra)

    result = CliRunner().invoke(pw2cw, ["ch4_sl"])

    assert result.exit_code == 1
    assert f"Error: {workdir / 'ch4_sl.cwin'}: " in result.output
    assert message in result.output


# ==================================================
def test_proj_min_too_large(make_case):
    make_case("ch4_sl", extra="proj_min = 1.1")

    result = CliRunner().invoke(pw2cw, ["ch4_sl"])

    assert result.exit_code == 1
    assert "proj_min = 1.1 is too large" in result.output
    assert "num_wann = 8" in result.output


# ==================================================
def test_postcw_without_mmn(make_case, monkeypatch):
    workdir = make_case("ch4_sl", extra="write_info_data = true")
    assert CliRunner().invoke(pw2cw, ["ch4_sl"]).exit_code == 0

    monkeypatch.chdir(workdir)
    result = CliRunner().invoke(postcw, ["ch4_sl"])

    assert result.exit_code == 1
    assert "requires ch4_sl.mmn" in result.output


# ==================================================
def test_exceptions_derive_from_builtin():
    assert issubclass(SymCWInputError, ValueError)
    assert issubclass(SymCWFileNotFoundError, FileNotFoundError)
    assert str(SymCWFileNotFoundError("eig", "a.eig")) == "cannot find the eig file: a.eig"
