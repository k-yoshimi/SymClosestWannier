"""
tests for reading input files compressed as seedname.ext.gz or seedname.ext.tar.gz.
"""

import io
import os
import gzip
import shutil
import tarfile

import numpy as np
import pytest
from scipy.io import FortranFile

from symclosestwannier.cw.cw_creator import cw_creator
from symclosestwannier.cw.mmn import Mmn
from symclosestwannier.cw.spn import Spn
from symclosestwannier.cw.uHu import UHu
from symclosestwannier.util.exceptions import SymCWFileNotFoundError, SymCWInputError
from symclosestwannier.util.input_file import open_input, input_path

from conftest import DATA_DIR, EXAMPLES_DIR, assert_hr_equal


# ==================================================
def compress(path, kind, member_dir=""):
    """
    replace path by path.gz or path.tar.gz.

    Args:
        path (pathlib.Path): file to compress.
        kind (str): "gz" or "tar.gz".
        member_dir (str, optional): directory of the member in tar archive.
    """
    if kind == "gz":
        with open(path, "rb") as src, gzip.open(f"{path}.gz", "wb") as dst:
            shutil.copyfileobj(src, dst)
    else:
        with tarfile.open(f"{path}.tar.gz", "w:gz") as tf:
            tf.add(path, arcname=os.path.join(member_dir, path.name))
    os.remove(path)


# ==================================================
@pytest.mark.parametrize("kind", ["gz", "tar.gz"])
def test_pw2cw_with_compressed_inputs(make_case, kind):
    seedname = "ch4_sl"
    workdir = make_case(seedname, extra="write_hr = true")
    for ext in ("eig", "amn", "nnkp"):
        compress(workdir / f"{seedname}.{ext}", kind, member_dir="data")

    cw_creator(seedname)

    assert_hr_equal(workdir / f"{seedname}_hr.dat.cw", os.path.join(DATA_DIR, f"{seedname}_hr.dat.cw"))


# ==================================================
def test_open_input_prefers_uncompressed(tmp_path):
    (tmp_path / "a.eig").write_text("plain\n")
    with gzip.open(tmp_path / "a.eig.gz", "wt") as fp:
        fp.write("gz\n")

    with open_input(str(tmp_path / "a.eig"), "eig") as fp:
        assert fp.read() == "plain\n"


# ==================================================
def test_open_input_single_member_with_other_name(tmp_path):
    (tmp_path / "other.txt").write_text("data\n")
    with tarfile.open(tmp_path / "a.eig.tar.gz", "w:gz") as tf:
        tf.add(tmp_path / "other.txt", arcname="other.txt")

    with open_input(str(tmp_path / "a.eig"), "eig") as fp:
        assert fp.read() == "data\n"


# ==================================================
def test_open_input_errors(tmp_path):
    with pytest.raises(SymCWFileNotFoundError, match="cannot find the eig file"):
        open_input(str(tmp_path / "a.eig"), "eig")

    for name in ("x.txt", "y.txt"):
        (tmp_path / name).write_text("")
    with tarfile.open(tmp_path / "a.eig.tar.gz", "w:gz") as tf:
        for name in ("x.txt", "y.txt"):
            tf.add(tmp_path / name, arcname=name)

    with pytest.raises(SymCWInputError, match="cannot find a.eig in"):
        open_input(str(tmp_path / "a.eig"), "eig")


# ==================================================
@pytest.mark.parametrize("kind", ["gz", "tar.gz"])
def test_input_path_removes_temporary_file(tmp_path, kind):
    (tmp_path / "a.spn").write_bytes(b"\x00\x01binary")
    compress(tmp_path / "a.spn", kind)

    with input_path(str(tmp_path / "a.spn"), "spn") as path:
        assert path != str(tmp_path / "a.spn")
        assert open(path, "rb").read() == b"\x00\x01binary"

    assert not os.path.exists(path)


# ==================================================
@pytest.mark.parametrize("kind", [None, "gz", "tar.gz"])
def test_mmn(tmp_path, kind):
    (tmp_path / "seed.mmn").write_text("header\n 1 1 1\n 1 1 0 0 0\n 0.5 0.25\n")
    if kind:
        compress(tmp_path / "seed.mmn", kind)

    mmn = Mmn(str(tmp_path), "seed", npar=1)

    np.testing.assert_allclose(mmn["Mkb"], [[[[0.5 + 0.25j]]]])


# ==================================================
def write_unformatted_spn(path):
    with FortranFile(str(path), "w") as f:
        f.write_record(np.frombuffer(b"header", dtype="S1"))
        f.write_record(np.array([1, 1], dtype=np.int32))
        f.write_record(np.array([1.0, 2.0, 3.0], dtype=np.complex128))


# ==================================================
@pytest.mark.parametrize("kind", [None, "gz", "tar.gz"])
@pytest.mark.parametrize("formatted", [True, False])
def test_spn(tmp_path, kind, formatted):
    if formatted:
        (tmp_path / "seed.spn").write_text("header\n 1 1\n 1.0 0.0\n 2.0 0.0\n 3.0 0.0\n")
    else:
        write_unformatted_spn(tmp_path / "seed.spn")
    if kind:
        compress(tmp_path / "seed.spn", kind)

    spn = Spn(str(tmp_path), "seed", formatted=formatted)

    np.testing.assert_allclose(np.array(spn["pauli_spn"]).reshape(3), [1.0, 2.0, 3.0])


# ==================================================
@pytest.mark.parametrize("kind", [None, "gz", "tar.gz"])
def test_uHu_formatted(tmp_path, kind):
    (tmp_path / "seed.uHu").write_text("header\n 1 1 1\n 0.5 -0.25\n")
    if kind:
        compress(tmp_path / "seed.uHu", kind)

    uHu = UHu(str(tmp_path), "seed", formatted=True)

    np.testing.assert_allclose(uHu["Hkb1b2"], [[[[[0.5 - 0.25j]]]]])


# ==================================================
def make_tar(path, members):
    """
    create tar.gz archive with members {name: text}; directories end with "/".
    """
    with tarfile.open(path, "w:gz") as tf:
        for name, text in members.items():
            info = tarfile.TarInfo(name.rstrip("/"))
            if name.endswith("/"):
                info.type = tarfile.DIRTYPE
                tf.addfile(info)
            else:
                data = text.encode()
                info.size = len(data)
                tf.addfile(info, io.BytesIO(data))


# ==================================================
def test_open_input_selects_member_by_name(tmp_path):
    make_tar(tmp_path / "a.eig.tar.gz", {"out/": "", "out/a.win": "win\n", "out/a.eig": "eig\n"})

    fp = open_input(str(tmp_path / "a.eig"), "eig")
    with fp:
        assert fp.read() == "eig\n"

    assert fp._tf.closed


# ==================================================
@pytest.mark.parametrize(
    "members",
    [
        {"x/a.eig": "1\n", "y/a.eig": "2\n"},  # ambiguous
        {"a.eig/": ""},  # no regular file
    ],
)
def test_open_input_rejects_archive(tmp_path, members):
    make_tar(tmp_path / "a.eig.tar.gz", members)

    with pytest.raises(SymCWInputError, match="cannot find a.eig in"):
        open_input(str(tmp_path / "a.eig"), "eig")


# ==================================================
def test_open_input_prefers_gz_over_tar(tmp_path):
    with gzip.open(tmp_path / "a.eig.gz", "wt") as fp:
        fp.write("gz\n")
    make_tar(tmp_path / "a.eig.tar.gz", {"a.eig": "tar\n"})

    with open_input(str(tmp_path / "a.eig"), "eig") as fp:
        assert fp.read() == "gz\n"


# ==================================================
@pytest.mark.parametrize("kind", ["gz", "tar.gz"])
@pytest.mark.parametrize("where", ["copy", "body"])
def test_input_path_cleanup_on_error(tmp_path, monkeypatch, kind, where):
    import symclosestwannier.util.input_file as input_file_module

    (tmp_path / "a.spn").write_bytes(b"data")
    compress(tmp_path / "a.spn", kind)

    created = []
    mkstemp = input_file_module.tempfile.mkstemp

    def record_mkstemp(*args, **kwargs):
        fd, path = mkstemp(*args, **kwargs)
        created.append(path)
        return fd, path

    monkeypatch.setattr(input_file_module.tempfile, "mkstemp", record_mkstemp)

    if where == "copy":

        def fail_copy(src, dst):
            raise OSError("disk full")

        monkeypatch.setattr(input_file_module.shutil, "copyfileobj", fail_copy)

    with pytest.raises(OSError, match="disk full"):
        with input_path(str(tmp_path / "a.spn"), "spn"):
            raise OSError("disk full")

    assert len(created) == 1
    assert not os.path.exists(created[0])


# ==================================================
def test_input_path_closes_descriptor_when_fdopen_fails(tmp_path, monkeypatch):
    import symclosestwannier.util.input_file as input_file_module

    (tmp_path / "a.spn").write_bytes(b"data")
    compress(tmp_path / "a.spn", "gz")

    closed = []
    monkeypatch.setattr(input_file_module.os, "fdopen", lambda fd, mode: (_ for _ in ()).throw(OSError("fdopen")))
    close = input_file_module.os.close
    monkeypatch.setattr(input_file_module.os, "close", lambda fd: (closed.append(fd), close(fd)))

    with pytest.raises(OSError, match="fdopen"):
        with input_path(str(tmp_path / "a.spn"), "spn"):
            pass

    assert len(closed) == 1
    assert sorted(p.name for p in tmp_path.iterdir()) == ["a.spn.gz"]


# ==================================================
def test_input_path_priority(tmp_path):
    (tmp_path / "a.spn").write_bytes(b"plain")
    with gzip.open(tmp_path / "a.spn.gz", "wb") as fp:
        fp.write(b"gz")
    make_tar(tmp_path / "a.spn.tar.gz", {"a.spn": "tar"})

    with input_path(str(tmp_path / "a.spn"), "spn") as path:
        assert open(path, "rb").read() == b"plain"

    os.remove(tmp_path / "a.spn")
    with input_path(str(tmp_path / "a.spn"), "spn") as path:
        assert open(path, "rb").read() == b"gz"


# ==================================================
def write_unformatted_uHu(path, num_bands=2, num_k=1, num_b=2):
    rng = np.random.default_rng(0)
    with FortranFile(str(path), "w") as f:
        f.write_record(np.frombuffer(b"header", dtype="S1"))
        f.write_record(np.array([num_bands, num_k, num_b], dtype=np.int32))
        for _ in range(num_k * num_b * num_b):
            f.write_record(rng.standard_normal(2 * num_bands * num_bands))


# ==================================================
@pytest.mark.parametrize("kind", ["gz", "tar.gz"])
def test_uHu_unformatted_compressed_equals_plain(tmp_path, kind):
    plain = tmp_path / "plain"
    packed = tmp_path / "packed"
    for d in (plain, packed):
        d.mkdir()
        write_unformatted_uHu(d / "seed.uHu")
    compress(packed / "seed.uHu", kind)

    expected = np.array(UHu(str(plain), "seed")["Hkb1b2"])
    result = np.array(UHu(str(packed), "seed")["Hkb1b2"])

    assert expected.shape == (1, 2, 2, 2, 2)
    np.testing.assert_array_equal(result, expected)


# ==================================================
def write_umat(path, num_wann=2, num_bands=3, num_k=1):
    rng = np.random.default_rng(1)

    def block(nrow, ncol):
        return "".join(f" {rng.standard_normal():.10f} {rng.standard_normal():.10f}\n" for _ in range(nrow * ncol))

    u = f"header\n {num_k} {num_wann} {num_wann}\n\n"
    u_dis = f"header\n {num_k} {num_wann} {num_bands}\n\n"
    for _ in range(num_k):
        u += " 0.0 0.0 0.0\n" + block(num_wann, num_wann) + "\n"
        u_dis += " 0.0 0.0 0.0\n" + block(num_bands, num_wann) + "\n"
    (path / "ch4_sl_u.mat").write_text(u)
    (path / "ch4_sl_u_dis.mat").write_text(u_dis)


# ==================================================
@pytest.mark.parametrize("kind", ["gz", "tar.gz"])
def test_umat_compressed_equals_plain(tmp_path, kind):
    from symclosestwannier.cw.umat import Umat

    plain = tmp_path / "plain"
    packed = tmp_path / "packed"
    for d in (plain, packed):
        d.mkdir()
        # seedname.win is read for dis_num_iter (= 0 for ch4_sl).
        shutil.copy(os.path.join(EXAMPLES_DIR, "ch4_sl", "ch4_sl.win"), d)
        write_umat(d)
    for name in ("ch4_sl_u.mat", "ch4_sl_u_dis.mat"):
        compress(packed / name, kind)

    expected = Umat(str(plain), "ch4_sl")
    result = Umat(str(packed), "ch4_sl")

    assert np.array(expected["Uk"]).shape == (1, 3, 2)
    for key in ("Uoptk", "Udisk", "Uk"):
        np.testing.assert_array_equal(np.array(result[key]), np.array(expected[key]))
