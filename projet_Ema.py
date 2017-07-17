# -*- coding: utf-8 -*-
"""
Created on Tue May 23 21:30:04 2017

@author: svincent
"""
from random import choice 
import binascii


def generate_random_key( n ):
    """
    This function generates a random encryption key. 
    The encryption key is a series of 0s and 1s. 
    The encryption key is of length n, where n is the size of the plaintext 
    message in binary format.
    For example if the message is kts, its binary format is 
    110101111101001110011. 
    The length of the binary message is 21. 
    From all possible encryption key of length 21, a valid encryption key is
    110011100011101100101000011000000011110100.
    Notice that if the legnth is 21, there exists 2^21=2097152 possible encryption keys
    """
    key  = []
    list = ['0','1']
    for i in range(n):
        key.append(choice(list))
    #return key, ''.join(key)
    print 'The encryption key for this session is:\n    ', ''.join(key)
    return ''.join(key)

#################################
#################################

def convert_message_to_bits( m ):
    """
    Given a plaintext message m, the function converts m into a binary format. 
    That is, m contains 0s and 1s only.
    For example if the plaintext message is kts then its binary format is:
    110101111101001110011
    """
    m_binary = []
    print 'The cleartext message is\n    ', m
    print 'Now converting the cleartext message into binary format...'
    for i in m:
        i_ascii  = ord(i)
        i_binary = '0'+format(i_ascii, 'b')
        m_binary.append(i_binary)
        #print '     character',i,'into ascii format',i_ascii,'and into binary format', i_binary
        print '     character',i,'into binary format', i_binary
    binary_message = ''.join(m_binary)
    print 'The cleartext message has now been converted into binary format\n    ',binary_message
    #print 'The length of the binary message is\n    ',len( binary_message )
    return binary_message
    
#################################
#################################    
    
def encrypt_message( m, k ):
    """
    Given a binary message m and an encryption key k of characters 0 and 1,
    the function returns the encrypted message. 
    Notice that m and k should be of the same length.
    """
    m_encrypted= []
    foo = zip( m, k )
    for i, j in foo:
        #print 'i:',i,'j:',j,'i xor j:',int(i,2) ^ int(j,2)
        m_encrypted.append(str(int(i,2) ^ int(j,2)))
    cipher = ''.join(m_encrypted) 
    print 'After xor-ing the encryption key to the binary message, our encrypted message is:\n    ', cipher
    return cipher

#################################
#################################

def get_text_from_binary( b ):
    """
    Given a binary message b, the function converts b into a text message.
    """
    foo = int(b, 2)
    return binascii.unhexlify('%x' % foo)

#################################
#################################

if __name__ == '__main__':
    #message = 'Bonjour ema'
    message = 'kts'
    binary_message    = convert_message_to_bits( message )
    encryption_key    = generate_random_key( len( binary_message ) )
    encrypted_message = encrypt_message(binary_message, encryption_key )
    #clear_message     = get_text_from_binary( binary_message )
    clear_message     = get_text_from_binary( encrypted_message )
    print '... that is:\n    ',clear_message
    