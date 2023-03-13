import json
from socket import socket
from math import ceil
from typing import Tuple

receive_byte_limit:int = 100000000#40960000
receive_timeout:int = 15
send_timeout:int = 15

def receive_object(socket:socket, timeout:bool=False) -> Tuple[object, bool]: 
    """Waits to receive an object on a socket with an optional timeout flag.
    Returns the received object if successful and a boolean value indicating
    success or failure.
    
    Args:
        socket:
        timeout:
    
    Returns:
        Tuple[object, bool] of (reconstructed_object, True) if successful else (None, False)
    """
    try:
        if timeout:
            socket.settimeout(receive_timeout)
        received_bytes = socket.recv(receive_byte_limit)
        if received_bytes == b'':
            # case of client disconnection,
            # return nothing and signal the
            # disconnect to the caller
            return (None, False)
        byte_count_bytes = received_bytes[:8]
        # this is the length of the incoming object json string
        # in bytes
        byte_count = int.from_bytes(byte_count_bytes, 'big', signed=False)

        # loop receiving bytes until the whole json string byte sequence has been received
        json_bytes = bytearray(received_bytes[8:])
        while len(json_bytes) < byte_count:
            next_bytes = socket.recv(receive_byte_limit)
            if next_bytes == b'':
                # case of client disconnection,
                # return nothing and signal the
                # disconnect to the caller
                return (None, False)
            json_bytes.extend(next_bytes)
        # convert bytearray to bytes
        json_bytes = bytes(json_bytes)

        json_string = json_bytes.decode("utf-8")
        #print("Received string:", json_string)
        # reconstruct the object sent from the client
        reconstructed_object = json.loads(json_string)
        # return object and signal successful receive
        return (reconstructed_object, True)
    except OSError: # except socket-related errors
        print("receive_object socket error")
        socket.close()
        return (None, False)
    finally:
        socket.settimeout(None)

def send_object(object:object, socket:socket, timeout:bool=False) -> bool:
    """ 
    add relevant docstring here 
    
    Args:
        object:
        socket:

    Returns:
        bool indicated send success or failure
    """
    try:
        if timeout:
            socket.settimeout(send_timeout)
        json_string = json.dumps(object)
        #print("Sending string:", json_string)
        json_bytes = bytes(json_string, encoding="utf-8")
        # encode length of bytes as an 8 byte (64 bit) number
        byte_count_bytes = len(json_bytes).to_bytes(8, 'big', signed=False)
        # send number of bytes in the json string bytes, encoded as a sequence of 8 bytes
        socket.sendall(byte_count_bytes)
        # send the bytes
        socket.sendall(json_bytes)
        return True
    except OSError as e:
        print("send_object socket error")
        print(e)
        socket.close()
        return False
    except TypeError as e:
        print("send_object type error")
        print(e)
        #print("Here here - TypeError")
        socket.close()
        return False
    finally:
        socket.settimeout(None)

#def send_object_with_new_connection(object, ip, port):
#    client_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
#    try:
#        # start server connection
#        client_socket.connect((ip, port))
#        send_object(object, client_socket)
#    finally:
#        client_socket.close()

#def receive_object_with_new_connection(ip, port):
#    client_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
#    try:
#        # start server connection
#        client_socket.connect((ip, port))
#        object, successp = receive_object(object, client_socket)
#        return object, successp
#    finally:
#        client_socket.close()
