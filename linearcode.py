import random
from fileutils import print_files


class LinearCode(object):

    def __init__(self, n, k, d, channel_error_probability=0.01):
        self.channel_error_probability = channel_error_probability
        self.n = n
        self.k = k
        self.d = d
        self.r = n - k
        self.h_matrix, self.g_matrix, self.a_matrix_transposed = get_source_matrices(
            n=self.n,
            r=self.r,
            d=self.d)
        self.syndrome_table = generate_syndrome_decoding_table(
            g_matrix=self.g_matrix,
            h_matrix=self.h_matrix,
            n=self.n,
            k=self.k,
            d=self.d)
        print_files(
            g_matrix=self.g_matrix,
            h_matrix=self.h_matrix,
            syndrome_table=self.syndrome_table,
            n=n,
            k=k)


def generate_syndrome_decoding_table(g_matrix, h_matrix, n, k, d):
    """
    Generates a syndrome decoding table,
    which maps a coset to its syndrome.

    :param g_matrix: a generator matrix.
    :param h_matrix: a parity check
    matrix.
    :param n: the length of the code word.
    :param k: a number of significant bits.
    :param d: a code distance.
    :return: syndrome_decoding_table -
    vectors have the same syndrome to the syndrome.
    """
    messages = 2 ** k
    syndrome_decoding_table = {0: [(
        multiply_matrices(
            matrix1=[i],
            columns_count1=k,
            matrix2=g_matrix,
            columns_count2=n)[0]) for i in range(messages)]}
    number = 2 ** n
    weights = get_binary_weight(number=number)
    j = 1
    cosets = 2 ** (n - k)
    for i in weights:
        matrix2 = transpose_matrix(matrix=[i], columns_count=n)
        syndrome_transposed = multiply_matrices(matrix1=h_matrix, columns_count1=n, matrix2=matrix2, columns_count2=1)
        syndrome = transpose_matrix(matrix=syndrome_transposed, columns_count=1)[0]
        if syndrome != 0:
            if get_hamming_weight(weights[j]) <= (d - 1) / 2:
                if syndrome not in syndrome_decoding_table:
                    syndrome_decoding_table[syndrome] = []
                syndrome_decoding_table[syndrome].append(weights[j])
            if j >= cosets:
                break
            j += 1
    return syndrome_decoding_table


def generate_i_matrix(size):
    return [(1 << (size - i - 1)) for i in range(size)]


def get_source_matrices(r, d, k):
    a_matrix = [0] * r
    xor_value = 0
    a_matrix_transposed = []
    for i in range(k):
        number = random.randrange(2 ** r)
        temp = xor_value ^ number
        if temp == 0: continue
        xor_value = temp
        fill_matrix_columns(
            matrix=a_matrix,
            index=i,
            number=number,
            n=k)
        a_matrix_transposed.append(number)
    return [add_identity_matrix(matrix1=generate_i_matrix(r), matrix2=a_matrix),
            add_identity_matrix(
                matrix1=a_matrix_transposed,
                matrix2=generate_i_matrix(k),
                lshift=r),
            a_matrix_transposed]


def transpose_matrix(matrix, columns_count):
    transposed_matrix = [0] * columns_count
    matrix_len = len(matrix)
    for i in range(matrix_len):
        num = matrix[i]
        fill_matrix_columns(
            matrix=transposed_matrix,
            index=i,
            number=num,
            n=matrix_len)
    return transposed_matrix


def add_identity_matrix(matrix1, matrix2, lshift=0):
    if len(matrix1) == len(matrix2):
        matrix_len = len(matrix1)
        matrix1 = [(matrix1[i] | matrix2[i] << lshift) for i in range(matrix_len)]
        return matrix1

    raise ValueError('The given matrices are different: len(identity_matrix): {0}, len(matrix): {1}'
                     .format(len(matrix1), len(matrix2)))


def fill_matrix_columns(matrix, index, number, n):
    r = len(matrix)
    matrix = [(matrix[j] | (1 << (n - 1 - index)) if number & (1 << (r - 1 - j)) else 0) for j in range(r)]


def get_binary_weight(number):
    length = len(bin(number)) - 2
    weights = {}
    for i in range(length):
        weights[i] = []
    for i in range(number):
        weight = get_hamming_weight(num=i)
        weights[weight].append(i)
    result = []
    for list_of_numbers in weights.values():
        for i in list_of_numbers:
            result.append(i)
    return result


def get_hamming_weight(num):
    length = len(bin(num)) - 2
    weight = 0
    for i in range(length):
        if num & 1 == 1:
            weight += 1
        num >>= 1
    return weight


def multiply_matrices(matrix1, columns_count1, matrix2, columns_count2):
    """
    Multiplies two matrices with the given number of columns.
    :return: the multuply result of two matrices.
    :raises: ValueError: if the count of
    columns of the first matrix doesn't equal to
    the number of lines of the second matrix.
    """
    m1_length = len(matrix1)
    m2_length = len(matrix2)

    if columns_count1 != len(matrix2):
        raise ValueError('The matrices are incompatible')

    transposed_matrix2 = transpose_matrix(matrix=matrix2, columns_count=columns_count2)
    result = [0] * m1_length

    for i in range(m1_length):
        for j in range(columns_count2):
            result[i] = result[i] | compute_parity(
                number=matrix1[i] & transposed_matrix2[j],
                length=m2_length) << (columns_count2 - j - 1)
    return result


def compute_parity(number, length):
    """
    Computes the parity, which is a number of 1s in
    the binary representation of a given number.
    Can be used to calculate the sum of all 1s.
    :return: 0 if the number of 1s is even, else 1.
    """
    number = number & (2 ** length - 1)
    i = 1
    while i <= length:
        number = number ^ (number >> i)
        i *= 2
    return number & 1
