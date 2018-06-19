# Connecting to the Vital-IT cluster

By now, you should have received a username and password to access the high performance computing cluster of the SIB (Vital-IT). You should use this login even though you already have a Vital-IT account, since we requested some privileges for the participants comfort:
 - 8 GB RAM by default
 - 24 h run-time
 - 12 dedicated nodes

In order to connect to the cluster and set up your working directory follow the next steps, depending on your OS:

## For Linux / Mac OSX users

In a terminal, type:

`ssh <username>@prd.vital-it.ch`

You will be prompted to input your user password:

- type in your password
- press "Enter"

You are in!

## For Windows users

You should first install a ssh client (e.g., [PuTTY](https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html)).

If you already have one, we assume you know how to use it.

If you do not have a ssh client, using PuTTY:
- download PuTTY for your platform (the 1st link is usually the best)
- set up PuTTY on your local machine ("Next")
- when prompted for PuTTY configuration, on the Connextion>SSH tab enter your credentials (e.g. `<username>@prd.vital-it.ch`), then hit "Open"
- when prompted for PuTTY security alert, trust the host to allow connection ("Yes")

You are in!

# Navigation

Go to [Set up & use working directory in Vital-IT](setup_vitalit.md)

Return to [Main page](README.md)
