import multiprocessing as mp
import tempfile
from pathlib import Path
from shutil import rmtree
from typing import Dict, Optional

from loguru import logger
from pyarrow.dataset import ParquetFileFragment

from . import pyarrow_reader as pr
from .helpers import set_cores
from .output import output_writer, tabix_accumulator


class ToTabix:
    """
    Convert Bed6 files to tabix formatted, indexed bed files
    """

    def __init__(
        self, bed6: str, out: str, cores: int, tmp_dir: Optional[str] = None
    ) -> None:
        self.bed6 = bed6
        self.out = out
        if tmp_dir is None:
            self._tmp = Path(self.out).parent / next(tempfile._get_candidate_names())
        else:
            if not Path(tmp_dir).exists():
                raise FileNotFoundError(f"Directory {tmp_dir} does not exists!")
            self._tmp: Path = Path(tmp_dir) / next(tempfile._get_candidate_names())
        logger.info(f"Using temporary directory: {self._tmp}")
        self._tmp.mkdir(parents=True, exist_ok=True)
        # set cores
        self.cores: int = set_cores(cores)
        self._pa_cores: int = max(1, int(self.cores / 2))
        self._rest_cores = max(1, self.cores - self._pa_cores)
        logger.debug(f"Using {self._pa_cores} for arrow and {self._rest_cores} for mp")

    def __enter__(self) -> "ToTabix":
        return self

    def __exit__(self, except_type, except_val, except_traceback):
        rmtree(self._tmp)  # clean up temp dir
        if except_type:
            logger.exception(except_val)

    def convert(self) -> None:
        """convert
        Convert input BED file to bgzipped, indexed tabix file
        """
        sorted_dict: Dict[str, str] = {}
        with pr.PartionedParquetReader(
            file_name=self.bed6,
            fformat="bed6",
            temp_dir=str(self._tmp),
            cores=self._pa_cores,
        ) as ppq:
            fragments = ppq.get_partitioned_fragments()
            with mp.Pool(processes=self._rest_cores) as pool:
                for chrom, fragment in fragments.items():
                    temp_file = str(
                        self._tmp / f"{next(tempfile._get_candidate_names())}.bed"
                    )
                    sorted_dict[chrom] = temp_file
                    pool.apply_async(bed_writer, args=(fragment, temp_file, chrom))
                # bed_writer(fragment=fragment, out=temp_file, chrom=chrom)
                pool.close()
                pool.join()
        tabix_accumulator(sorted_dict, str(self._tmp), self.out, "bed")


def bed_writer(fragment: ParquetFileFragment, out: str, chrom: str) -> None:
    """bed_writer
    Worker function, sort each parquet fragment by chromosome start and write to a temp. file
    Args:
        fragment: ParquetFileFragment, parquet file fragment for chromosome
        out: str, path to write bed file
        chrom: str, chromosome name
    Returns:
        None
    """
    logger.debug(
        f"Process id: {mp.current_process().pid}, chromosome: {chrom}, temp. file: {out}, max. chunksize: {10000:,}"
    )
    bwriter = output_writer(out, use_tabix=False, preset="bed")
    with bwriter(out) as _ow:
        sorted_batch = (
            fragment.to_table().sort_by("chromStart").to_batches(max_chunksize=10000)
        )
        logger.debug(f"chromosome: {chrom}, chunks: {len(sorted_batch)}")
        for batch in sorted_batch:
            for entry in batch.to_pylist():
                _ow.write(
                    f"{chrom}\t{entry['chromStart']}\t{entry['chromEnd']}\t{entry['name']}\t{entry['score']}\t{entry['strand']}\n"
                )
