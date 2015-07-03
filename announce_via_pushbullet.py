def main():

    from yapbl import PushBullet
    import os

    assert(os.path.isfile('access_token.txt'))
    with open('access_token.txt') as f:
        access_token = f.readline().rstrip()

    print(access_token)

    p = PushBullet(access_token)
    p.push_note('the simulation is finished',
                os.getcwd())

if __name__ == '__main__':

    main()

    
