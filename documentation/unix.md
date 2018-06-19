# You'll be logging to a UNIX-based OS

The server you are using as working platform is based on [CentOS](https://en.wikipedia.org/wiki/CentOS): you can use [Unix](https://en.wikipedia.org/wiki/Unix) to communicate with it.

We assume that you have some sort of familiarity with Unix, since you had the opportunity to take the [SIBWeb Unix Fundamentals](https://edu.isb-sib.ch/course/view.php?id=82) e-learning module.

In case you would like a refresh, we provide you a [Unix cheat sheet](http://cheatsheetworld.com/programming/unix-linux-cheat-sheet/): the section *File system* should be enough to satisfy the purposes of this tutorial.

We try to be as clear as possible, but in case we use technical concepts you don't fully understand, we provide you a [Unix vocabulary](https://www.cs.bu.edu/teaching/unix/reference/vocab.html). You can of course also ask us directly (we don't bite ☺ ).

# Some basic commands

Once you're logged to the server, you should find yourself at home.

To check your location, you can use the command `pwd`: this will return the absolute path to your current folder.

```
$ pwd
/home/<username>
```

To create a new directory, you can use the command `mkdir`: this will create a new folder from within your current folder, unless otherwise specified.

```
mkdir <foldername>            # will create <foldername> with absolute path: /home/<username>/<foldername>
mkdir /path/to/<foldername>   # will create <foldername> with absolute path: /path/to/<foldername>
```

To move to a directory, you can use the command `cd`: this will move to a new folder from within your current folder, unless otherwise specified. You can thus use `cd` command to browse through directories, just like you would do with your favorite file manager.

```
cd <foldername>               # will move to <foldername> with absolute path: /home/<username>/<foldername>
cd /path/to/<foldername>      # will move to <foldername> with absolute path: /path/to/<foldername>
cd ~                          # will move you to your home directory, equivalent of: cd /home/<username>
```

To remove a file or a directory, you can use the command `rm`. Please be aware that once you remove a file/directory, this one cannot be recovered anymore (it's not like a bin in GUI).

```
rm file                       # will remove the selected file
rm -r folder/                 # will remove the selected folder
```

To rename a file or a directory, you can use the command `mv`. This command can also be used to move a file through directories.

```
mv file new_file              # will rename file to new_file
mv file ../                   # will move file to the higher level directory
```

To list the files in a given directory, you can use the command `ls`: this will print a simple list of files and folder names present in the directory from which you execute the command. For a more detailed list, also including e.g., size and permissions, you might want to try `ls -lh`.

```
cd ~
ls -lh
```

To download data from the internet given the location on a web server, you can use the command `wget`: this will allow you to retrieve the content via HTTP(/S) and FTP. You can follow the download status checking the standard output flow (aka, your terminal).

```
wget <url>
```

To copy some files from a remote host (e.g., the Amazon machine) to the local host (e.g., your PC), you can use the the command `scp` - it uses ssh for data transfer and provides the same authentication as `ssh`.

```
scp username@remotehost:/path/to/file.txt /path/to/local/folder
```

# A very cool feature: manual pages

Manual pages (aka, man pages) are useful to access the documentation for a program you're interested in, ideally covering all the features that the program offers. For example, try it on one of the commands you just used. Note that if the command you are looking for is a bash builtin, `man` will redirect to the builtin man page:

```
man ls
man cd
```

By default, this will open a terminal pager program (e.g., `less`), so:
* to scroll the manual you just use the up/down arrows
* to quit it you need to hit the letter `q`
* to perform a full-text search you type your keyword with `/<keyword>`
* to navigate to the next highlighted instance with the letter `n`


Don't hesitate to learn more about Unix commands using man pages!
