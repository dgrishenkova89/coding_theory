import unittest
import yaml
import random

from fileutils import read_file_to_dict, read_file_to_list
from correctingcode import CorrectingCode, \
    is_low_bound, \
    calculate_parity_check_error, \
    encode, decode
from linearcode import get_source_matrices, \
    get_hamming_weight, \
    compute_parity, \
    generate_syndrome_decoding_table, \
    transpose_matrix, \
    multiply_matrices


class LinearCodeTest(unittest.TestCase):
    n = 21
    k = 6
    r = 15
    t = 3
    d = 2 * t + 1

    def test_init_l_c_timeout(self):
        with self.assertRaises(TimeoutError):
            CorrectingCode(n=self.n, k=self.k, d=self.d, channel_error_probability=0.01)

    def test_init_e_c_c(self):
        CorrectingCode(n=self.n, k=self.k, d=self.d, channel_error_probability=0.01)

    def test_is_gilbert_varshamov_bound_must_assert(self):
        assert not is_low_bound(n=10, k=5, d=4)
        assert is_low_bound(n=12, k=4, d=4)

    def test_calculate_parity_check_error(self):
        n = 10
        t = 1
        p = 0.01
        assert round(calculate_parity_check_error(n=n, t=t, p=p), 5) == 0.00427

    def test_generate_a_matrix_and_a_matrix_transposed(self):
        h_matrix, g_matrix, a_matrix_transposed = get_source_matrices(
            r=self.r,
            d=self.d,
            k=self.k)
        print("\r\nH matrix: ", h_matrix)
        print("G matrix: ", g_matrix)
        print("Transposed matrix: ", a_matrix_transposed)
        linear_row = [(row != 0) for row in h_matrix]
        bound_hamming_weight = [(get_hamming_weight(i) >= self.d - 1) for i in a_matrix_transposed]

        assert len(linear_row) != 0
        assert len(bound_hamming_weight) != 0

    def test_compute_parity(self):
        assert compute_parity(number=0b11111111111, length=11) == 1
        assert compute_parity(number=0b111111111111, length=3) == 1

    def test_generate_syndrome_decoding_table(self):
        n = 6
        r = 3
        k = n - r
        d = 1
        h_matrix, g_matrix, a_matrix_transposed = get_source_matrices(r=r, d=d, k=k)
        assert len(h_matrix) != 0
        assert len(g_matrix) != 0
        assert len(a_matrix_transposed) != 0
        syndrome_decoding_table = generate_syndrome_decoding_table(
            g_matrix=g_matrix,
            h_matrix=h_matrix,
            n=n,
            k=k,
            d=d)
        for key in syndrome_decoding_table.keys():
            if key != 0:
                for word in syndrome_decoding_table[key]:
                    assert get_hamming_weight(word) <= self.t
                    assert word != 0

    def test_read_file_to_dict(self):
        config = yaml.safe_load(open('config-test.yml'))
        dictionary = read_file_to_dict(name=config['syndrome-decoding-test'])
        n = 10
        k = 3
        r = n - k
        zero_word_counter = 0
        word_counter = 0
        for key in dictionary.keys():
            for word in dictionary[key]:
                if word == 0:
                    zero_word_counter += 1
                word_counter += 1
        assert zero_word_counter <= 1
        assert word_counter > 2 ** r

    def test_code(self):
        config = yaml.safe_load(open('config-test.yml'))
        for i in range(2 ** self.k):
            code, distorted_code, error = encode(
                coder_file=config['coder-generator-test'],
                message=i,
                m_length=self.k,
                error=74)  # 000000000000001001010
            assert transpose_matrix(
                matrix=multiply_matrices(
                    matrix1=read_file_to_list(config['decoder-parity-check-test']),
                    columns_count1=self.n,
                    matrix2=transpose_matrix([code], self.n),
                    columns_count2=1
                ),
                columns_count=1)[0] == 0
        for i in range(2 ** self.k):
            code, distorted_code, error = encode(
                coder_file=config['coder-generator-test'],
                message=i,
                m_length=self.k)
            if get_hamming_weight(error) <= self.t:
                assert transpose_matrix(
                    matrix=multiply_matrices(
                        matrix1=read_file_to_list(config['decoder-parity-check-test']),
                        columns_count1=self.n,
                        matrix2=transpose_matrix([code], self.n),
                        columns_count2=1
                    ),
                    columns_count=1)[0] == 0

    def test_code_throw_exception(self):
        config = yaml.safe_load(open('config-test.yml'))
        with self.assertRaises(ValueError):
            encode(
                coder_file=config['coder-generator-test'],
                message=1,
                m_length=2,
                error=74)

    def test_decode(self):
        message = 56
        config = yaml.safe_load(open('config-test.yml'))
        decoded_message = decode(
            parity_check_file=config['decoder-parity-check-test'],
            n=self.n,
            syndrome_file=config['syndrome-decoding-test'],
            distorted_code=int('110010000110111000100', 2)
        )
        assert decoded_message == message

    def test_code_decode(self):
        config = yaml.safe_load(open('config-test.yml'))
        number_of_unfixed_errors = 0
        number_of_repetitions = 100
        for i in range(number_of_repetitions):
            message = random.randrange(2 ** self.k)
            rand_error = get_hamming_weight(self.n, int((self.d - 1) / 2))
            code, distorted_code, error = encode(
                coder_file=config['coder-generator-test'],
                message=message,
                m_length=self.k,
                error=rand_error)
            assert rand_error == error
            decoded = decode(
                parity_check_file=config['decoder-parity-check-test'],
                n=self.n,
                syndrome_file=config['syndrome-decoding-test'],
                distorted_code=distorted_code
            )
            if decoded != message:
                number_of_unfixed_errors += 1
        assert number_of_unfixed_errors == 0