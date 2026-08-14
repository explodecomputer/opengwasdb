"""Dense layout defaults."""

DEFAULT_CHUNK_SHAPE = (1000, 1000)
DEFAULT_COMPRESSOR = {
    "library": "numcodecs.Blosc",
    "cname": "zstd",
    "clevel": 3,
    "shuffle": "bitshuffle",
}
DEFAULT_DTYPE = "float16"
