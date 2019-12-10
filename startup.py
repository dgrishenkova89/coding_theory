import yaml
from correctingcode import CorrectingCode, get_random_hamming_weight, encode, decode

from fileutils import remove_files


def startup():
    remove_files()
    r, n, t = [int(x) for x in
               input("Enter number of check bits, length of a code word, number of error to correct: ").split(',')]
    k = n - r
    d = 2 * t + 1
    CorrectingCode(n=n, k=k, d=d)
    message = int(input("Enter a message: "), 2)
    error_value = int(input("Enter an error in binary, 0 if you want to skip: "), 2)
    if error_value == 0:
        error_value = get_random_hamming_weight(n, t)
    config = yaml.safe_load(open('config.yml'))
    codeword, distorted_codeword, error_value = encode(
        coder_file=config['Code_G_Matrix'],
        message=message,
        m_length=k,
        error=error_value)
    print("Code word: {0:>0{a}b}".format(codeword, a=n))
    print("Distorted code word: {0:>0{a}b}".format(distorted_codeword, a=n))
    print("Error vector: {0:>0{a}b}".format(error_value, a=n))
    distorted_codeword = int(input("Enter a distorted codeword: "), 2)
    decoded_message = decode(
        parity_check_file=config['Decoder_H_Matrix'],
        n=n,
        syndrome_file=config['Syndrome_Decoding'],
        distorted_code=distorted_codeword)
    print("Decoded message: {0:>0{a}b}".format(decoded_message, a=n - r))
