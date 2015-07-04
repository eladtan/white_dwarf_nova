def main():

    from yapbl import PushBullet
    import os
    import glob

    assert(os.path.isfile('access_token.txt'))
    with open('access_token.txt') as f:
        access_token = f.readline().rstrip()

    p = PushBullet(access_token)

    mp4_list = glob.glob('*.mp4')
    if len(mp4_list)>0:
        for fname in mp4_list:
            p.push_file(open(fname,'rb'))
    p.push_note('the simulation is finished',
                os.getcwd())

if __name__ == '__main__':

    main()

    
