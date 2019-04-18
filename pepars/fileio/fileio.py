import urllib.request
import shutil
import os


def download_remote_file(file_URL, file_path, skip_if_exists=True):
    """
    Download a file from the given URL to the given local file_path,
    if it doesn't already exist

    :param file_URL: An http or https URL to a direct file download
    :param file_path: The path of where to download the file to locally
    :param skip_if_exists: Whether to skip if this file already exists

    :return: None
    """

    if os.path.exists(file_path):
        if skip_if_exists:
            return
        else:
            os.remove(file_path)

    tmp_file_path = file_path + ".tmp"

    if os.path.exists(tmp_file_path):
        os.remove(tmp_file_path)

    directory_path = os.path.dirname(file_path)

    if directory_path:
        if not os.path.exists(directory_path):
            os.makedirs(directory_path)

    # Download the file from `url` and save it locally under `file_name`:
    with urllib.request.urlopen(file_URL) as response,\
            open(tmp_file_path, "wb") as out_file:
        shutil.copyfileobj(response, out_file)

    shutil.move(tmp_file_path, file_path)


