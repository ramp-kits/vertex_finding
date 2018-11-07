from __future__ import print_function

import os
import tarfile

try:
    from urllib.request import urlretrieve
except ImportError:
    from urllib import urlretrieve

#helper function to unpack the data
def unpack_tar(filename, type):
    tar = tarfile.open(filename)
    for member in tar.getmembers():
      #check if member is file
      if member.isreg():
        name = member.name
        if name.endswith('.json'):
          #remove prefix
          name = name.split('/')[1]
          member.name = name
          tar.extract(member, './data/' + type)
    tar.close()


URLBASE = 'https://storage.ramp.studio/vertex_finding/{}'
DATA = [
    '28kEvents.tar.gz', '30kEvents.tar.gz']



def main(output_dir='data'):
    filenames = DATA 
    urls = [URLBASE.format(filename) for filename in filenames]

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    if not os.path.exists(output_dir +'/train'):
        os.mkdir(output_dir +'/train')
    if not os.path.exists(output_dir +'/test'):
        os.mkdir(output_dir +'/test')

    # notfound = []
    for url, filename in zip(urls, filenames):
        output_file = os.path.join(output_dir, filename)

        if os.path.exists(output_file):
            continue

        print("Downloading from {} ...".format(url))
        urlretrieve(url, filename=output_file)
        print("=> File saved as {}".format(output_file))
    print("Unpacking tar files. This may take a while, please be patient.")
    unpack_tar('data/30kEvents.tar.gz', 'train')
    unpack_tar('data/28kEvents.tar.gz', 'test')


if __name__ == '__main__':
    test = os.getenv('RAMP_TEST_MODE', 0)

    if test:
        print("Testing mode, not downloading any data.")
    else:
        main()
