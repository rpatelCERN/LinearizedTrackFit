__author__ = 'demattia'


def combination_index(layers, radius):
    """
    Evaluates the combination index from a list of layers/disks and the radii (usef only for the disks)
    """
    bit_list = [0]*32
    for l, r in zip(layers, radius):
        bit_list[l] = 1
        if l > 10 and r > 61:
            bit_list[l+10] = 1
    combination_index_value = 0
    # Convert the list of bits to integer
    for bit in reversed(bit_list):
        combination_index_value = (combination_index_value << 1) | bit
    return combination_index_value


class CombinationsGenerator:
    def __init__(self):
        self.combs = {}
        self.comb_length = 0
        self.combinations_index = 0
        self.combination_ = []
        # # Generate all combinations for 6/8 and for 5/8 when there are 8 and store them.
        # indexes = range(8)
        # self.combinations(indexes, 7)
        # self.combinations(indexes, 6)
        # self.combinations(indexes, 5)
        #
        # # Generate all combinations for 6/7 and for 5/7 when there are 7 and store them.
        # indexes = range(7)
        # self.combinations(indexes, 6)
        # self.combinations(indexes, 5)

        # Generate all default combinations for 5/6 when there are 6 and store them.
        indexes = range(0, 6)
        self.combinations(indexes, 5)

    def combinations_size(self, combinations_index):
        return len(self.combs[combinations_index])

    def combination(self, index, combinations_index):
        return self.combs[combinations_index][index]

    def combinations(self, indexes, comb_length):
        self.comb_length = comb_length
        # The index of this set of combinations is the length in the input vector.
        self.combinations_index = len(indexes)

        if self.combinations_index not in self.combs:
            self.combs[self.combinations_index] = []
        # if comb_length > self.combinations_index:
        #     except
        if comb_length == self.combinations_index:
            self.combs[self.combinations_index].append(indexes)
            return self.combs[self.combinations_index]
        self.combination_ = [0]*self.comb_length
        self.generate_combinations(indexes, 0, 0)
        return self.combs[self.combinations_index]

    def generate_combinations(self, indexes, start, position_in_combination):
        if position_in_combination == self.comb_length:
            # Use the list() to make a copy
            self.combs[self.combinations_index].append(list(self.combination_))
        else:
            for i in range(start, self.combinations_index - self.comb_length + position_in_combination + 1):
                self.combination_[position_in_combination] = indexes[i]
                self.generate_combinations(indexes, i+1, position_in_combination+1)
