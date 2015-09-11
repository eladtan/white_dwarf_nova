def main():
    """
    Requires youtube uploader script from
    https://github.com/tokland/youtube-upload
    """

    from subprocess import Popen, PIPE
    import glob
    from sh import git

    yt_ids = []
    for fname in glob.glob('*.webm'):
        title = fname.replace('.webm','').replace('_',' ')
        command = 'youtube-upload --title="'+title+'" '+fname
        p = Popen(command,stdout=PIPE,shell=True)
        out = p.communicate()
        yt_ids.append(str(out[0].rstrip()).replace("b'",'').replace("'",''))
    readme_content = '# White dwarf nova\n'
    for idd in yt_ids:
        readme_content += '[![IMAGE ALT TEXT HERE](http://img.youtube.com/vi/'+idd+'/0.jpg)](http://www.youtube.com/watch?v='+idd+')\n'
    with open('README.md','w') as f:
        f.write(readme_content)
    git.add('README.md')
    git.commit(m='update videos')
    git.push()

if __name__ == '__main__':

    main()
