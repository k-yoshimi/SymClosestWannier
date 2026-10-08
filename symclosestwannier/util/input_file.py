"""
open input files given as file_name, file_name.gz or file_name.tar.gz.
"""

import os
import io
import gzip
import shutil
import tarfile
import tempfile
import contextlib

from symclosestwannier.util.exceptions import SymCWFileNotFoundError, SymCWInputError


# ==================================================
def _tar_member(tf, file_name, archive):
    """
    member of tar archive for file_name: the file with the same base name, or the only file in the archive.

    Args:
        tf (tarfile.TarFile): tar archive.
        file_name (str): file name without .tar.gz.
        archive (str): file name of the archive.

    Returns:
        tarfile.TarInfo: member.
    """
    files = [m for m in tf.getmembers() if m.isfile()]
    same_name = [m for m in files if os.path.basename(m.name) == os.path.basename(file_name)]

    if len(same_name) == 1:
        return same_name[0]
    if len(files) == 1:
        return files[0]

    raise SymCWInputError(f"cannot find {os.path.basename(file_name)} in {archive}.")


# ==================================================
class _TarTextFile(io.TextIOWrapper):
    """
    text file of a tar member, which closes the archive when closed.
    """

    def __init__(self, tf, member):
        self._tf = tf
        super().__init__(tf.extractfile(member))

    def close(self):
        super().close()
        self._tf.close()


# ==================================================
def open_input(file_name, kind):
    """
    open file_name, file_name.gz or file_name.tar.gz as a text file (in this order of priority).

    Args:
        file_name (str): file name.
        kind (str): kind of file used in the error message, e.g. "eig".

    Returns:
        file object: text file object, to be closed by the caller.
    """
    if os.path.exists(file_name):
        return open(file_name, "r")

    if os.path.exists(file_name + ".gz"):
        return gzip.open(file_name + ".gz", "rt")

    if os.path.exists(file_name + ".tar.gz"):
        archive = file_name + ".tar.gz"
        tf = tarfile.open(archive, "r:gz")
        try:
            return _TarTextFile(tf, _tar_member(tf, file_name, archive))
        except BaseException:
            tf.close()
            raise

    raise SymCWFileNotFoundError(kind, file_name)


# ==================================================
@contextlib.contextmanager
def input_path(file_name, kind):
    """
    path of file_name, file_name.gz or file_name.tar.gz as an uncompressed file.
    compressed files are extracted to a temporary file, which is removed afterwards.
    used for libraries which need a file name, e.g. fortio for unformatted files.

    Args:
        file_name (str): file name.
        kind (str): kind of file used in the error message, e.g. "spn".

    Yields:
        str: path of the uncompressed file.
    """
    if os.path.exists(file_name):
        yield file_name
        return

    # resources are closed (and the temporary file is removed) in reverse order, also on errors.
    with contextlib.ExitStack() as stack:
        if os.path.exists(file_name + ".gz"):
            src = stack.enter_context(gzip.open(file_name + ".gz", "rb"))
        elif os.path.exists(file_name + ".tar.gz"):
            archive = file_name + ".tar.gz"
            tf = stack.enter_context(tarfile.open(archive, "r:gz"))
            src = stack.enter_context(tf.extractfile(_tar_member(tf, file_name, archive)))
        else:
            raise SymCWFileNotFoundError(kind, file_name)

        fd, tmp = tempfile.mkstemp(prefix=os.path.basename(file_name) + ".")
        stack.callback(os.remove, tmp)
        try:
            dst = os.fdopen(fd, "wb")
        except BaseException:
            os.close(fd)
            raise
        with dst:
            shutil.copyfileobj(src, dst)

        yield tmp
