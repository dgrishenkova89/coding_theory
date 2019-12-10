import random
import math

from fileutils import read_file_to_list, read_file_to_dict
from linearcode import LinearCode, transpose_matrix, multiply_matrices


class CorrectingCode(LinearCode):

    def __init__(self, n, k, d, channel_error_probability=0.01):
        """

        r: a number of check bits.
        n: a codeword length.
        t: a number of errors, which the code should correct.
        """
        self.channel_error_probability = channel_error_probability
        self.n = n
        self.k = k
        self.d = d
        self.r = n - k
        if not is_low_bound(n=self.n, k=self.k, d=self.d):
            raise ValueError(
                'The given n: {0}, k: {1} and d: {2} aren\'t compliant with the Gilbert-Varshamov bound.'.format(n, k, d))

        LinearCode.__init__(
            self=self,
            n=self.n,
            k=self.k,
            d=self.d,
            channel_error_probability=self.channel_error_probability)

        print('The decoder error probability is {0:.2%}'.format(
            calculate_parity_check_error(n=self.n, t=self.t, r=self.r, p=self.channel_error_probability)))


def is_low_bound(n, k, d):
    """
    Check the specified parameters aren't complaint with the
    Gilbert-Varshamov bound.
    """
    result = 1
    for i in range(1, d - 1):
        result += math.factorial(n - 1) / math.factorial(i) / math.factorial(n - 1 - i)
    is_in_bound = result < 2 ** (n - k)
    return is_in_bound


def calculate_parity_check_error(n, t, p):
    sum = 0
    for i in range(t + 1):
        r = n - i
        binomial_coefficient = math.factorial(n) / math.factorial(i) / math.factorial(r)
        sum += binomial_coefficient * p ** i * (1 - p) ** r
    probability = 1 - sum
    print("Probability: {0:.2%}".format(probability))
    return probability


def encode(coder_file, message, m_length, error=0):
    generator_matrix = read_file_to_list(name=coder_file)
    max_line = max(generator_matrix)
    number_of_columns_g_m = len(bin(max_line)) - 2
    gen_matrix_length = len(generator_matrix)
    if gen_matrix_length != m_length:
        raise ValueError('The given len(generator_matrix) == {0}'
                         ' and m_length == {1} are incompatible.'
                         ' {0} != {1}.'.
                         format(gen_matrix_length, m_length))
    if error == 0:
        error = random.randrange(2 ** number_of_columns_g_m)
    code = multiply_matrices(
        matrix1=[message],
        columns_count1=m_length,
        matrix2=generator_matrix,
        columns_count2=number_of_columns_g_m)[0]
    return code, code ^ error, error


def decode(parity_check_file, n, syndrome_file, distorted_code):
    parity_check_matrix = read_file_to_list(name=parity_check_file)
    syndrome_decoding_table = read_file_to_dict(name=syndrome_file)
    syndrome = transpose_matrix(
        matrix=multiply_matrices(
            matrix1=parity_check_matrix,
            columns_count1=n,
            matrix2=transpose_matrix(
                matrix=[distorted_code],
                columns_count=n),
            columns_count2=1),
        columns_count=1)[0]
    return syndrome_decoding_table[0].index(distorted_code ^ syndrome_decoding_table[syndrome][0])


def get_random_hamming_weight(length, weight):
    if weight > length:
        raise ValueError('The weight shouldn\'t be greater than the length: {0} > {1}'.format(weight, length))
    i = 0
    result = 0
    while True:
        if i == weight:
            return result
        shift = random.randrange(length)
        power_of_two = 1 << shift
        if power_of_two & result == power_of_two:
            continue
        result |= power_of_two
        i += 1