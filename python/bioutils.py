def partition_lst(lst, part):
    """Partition lst into part separate lists of items.
    :param lst: list to be partitioned
    :param part: the number of lists to generate"""

    partition_size = len(lst)/part
    """Create a list of file lists."""
    partition_list = list()
    index = 0

    """Partition the list by multiplying an index by the partition size to
    get the starting index, and add the partition_size to get the end
    index. If the number of files does not divide evenly by NUM_CORES,
    we must add the tail end of the file list outside of the loop."""
    while (index < part):
        start_index = index*partition_size
        partition_list.append(lst[start_index:start_index+partition_size])
        index += 1
    if (len(lst) % part):
        """Allocate the leftover files over all the processors."""
        start_index = index*partition_size
        leftover_list = lst[start_index:len(lst)]
        mod_index = 0
        for f in leftover_list:
            """Modulo the index to wrap around the file_list."""""
            partition_list[mod_index % len(lst)].append(f)
            mod_index += 1

    return partition_list

