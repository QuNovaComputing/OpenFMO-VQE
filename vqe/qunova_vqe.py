import copy
from http import client
import socket, ssl
import os

import numpy as np

from abc import ABC
from typing import Union, List, Tuple, Optional

#import os
from subprocess import call, Popen
from time import sleep

from object_transfer import send_object, receive_object

_port:int = 40675
_ip:str = "localhost"
_ssl_certificate_path:str = os.path.join(os.path.dirname(os.path.abspath(__file__)), "vqe_server_cert.pem")

# NOTE: unverified context can be changed to default context once the server has a 
# fixed hostname and a certificate authority can verify the business details on the certificate.
# Must ensure ssl options match server's ssl options
_ssl_context = ssl._create_unverified_context(cafile=_ssl_certificate_path,
                                             #certfile=client_certificate_path, #if client needs to send its own certificate, give it here
                                             protocol=ssl.PROTOCOL_TLS_CLIENT,
                                             cert_reqs=ssl.CERT_REQUIRED) #cert_reqs=ssl.CERT_REQUIRED forces server to authenticate itself
_client_socket:socket = None



def establish_pulsar_connection(username:str, api_key:str, persistent_session_p=False) -> None:
    """
    Establishs a connection with the pulsar VQE server.

    Args:
        username:string username
        api_key: string API key
        persistent_session_p: boolean. When True, tasks submitted will not be tied to the
                                       current connection session and will not be canceled
                                       when the client disconnects from the server.

    Returns:
        None
    """
    if type(username) != str:
        raise ValueError("Username must be a string.")
    elif type(api_key) != str:
        raise ValueError("API key must be a string.")

    global _client_socket
    if _client_socket is not None:
        _client_socket.close()
    attempt_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    attempt_socket = _ssl_context.wrap_socket(attempt_socket)
    # connect to the server
    attempt_socket.connect((_ip, _port))
    # send the username and api key to the server
    authentication_data = {'username':username, 'api_key':api_key, 'persistent_session_p':persistent_session_p}
    successp = send_object(authentication_data, attempt_socket, timeout=True) #timeout to avoid client waiting forever
    if not successp:
        raise Exception("Could not send API authentication details to server.")
    server_response, successp = receive_object(attempt_socket) # timeout not required, since server disconnect will also trigger a response
    if not successp:
        raise Exception("Did not receive validation response from server.")
    if server_response["result"] != "success":
        raise Exception(server_response["reason"])
    # finally, set the global connection socket variable
    # to the successfully connected and authenticated attempt_socket
    # object    
    _client_socket = attempt_socket


def _interact_with_server(to_send:object) -> object:
    """
    Allows sending an object to the server and then receiving one in return. Also handles ensuring the connection stays alive.

    Args:
    
    """
    if _client_socket is None:
        raise Exception(f"Please call establish_pulsar_connection from the {__name__} module with your username and api key before attempting to use functionality requiring server access.")
    
    #TODO: add exception handling in this function for when the server interaction fails in either direction
    # send 1 object to the server and receive one in return
    successp = send_object(to_send, _client_socket, timeout=True) # in cases where the client sends information, the client handler will always be ready to receive it, so timeout to prevent clients waiting forever
    assert successp

    # loop receiving and responding to keepalive messages and finally receive the returned result
    try:
        while True:
            # wait for server response
            returned_object, successp = receive_object(_client_socket)
            assert successp
            if 'message_type' in returned_object and returned_object['message_type'] == "keepalive":
                # if asked to do a keepalive response, send one
                send_object({'message_type':'keepalive'}, _client_socket, timeout=True)
            else:
                # otherwise we were just sent the result, so stop looping
                break
    except OSError as e:
        print("Socket error interacting with server.")
        print(e)
        raise e

    return returned_object


# Define generators of functions that make basic server queries
def _generate_id_query_function(operation):
    def query_function(task_id):
        returned_object = _interact_with_server({"operation":operation, "task_id":task_id})
        if returned_object["result"] == "failure":
            raise ValueError(returned_object["reason"])
        else:
            return returned_object["result"]
    return query_function


# use the first generator to define the basic query functions for specific tasks
get_result = _generate_id_query_function("get_result") #function blocks until task result is ready, then returns it
def _submit_task_and_get_result(task_data):
    returned_object = _interact_with_server(task_data)
    if returned_object["result"] == "failure":
        raise ValueError(returned_object["reason"])
    else:
        task_id = returned_object["task_id"]
        task_result = get_result(task_id)
        return task_result

def call_vqe(contents):

    contents['oei'] = contents['oei'].tolist()
    contents['s'] = contents['s'].tolist()
    contents['c'] = contents['c'].tolist()
    contents['tei'] = contents['tei'].tolist()
    contents['env_oei'] = contents['env_oei'].tolist()
    #contents['moons'] = contents['moons'].tolist()

    object = {'operation':'call_vqe',
              'contents':contents}

    result = _submit_task_and_get_result(object)
    
    new_res = (result[0], result[1], result[2], result[3], np.array(result[4]))
    
    return new_res
