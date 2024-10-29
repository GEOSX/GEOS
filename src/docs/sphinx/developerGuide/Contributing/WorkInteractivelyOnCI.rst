
.. _WorkingInteractivelyOnCI:

********************************************
How to work interactively on the CI Machines
********************************************

When developing with GEOS, developers may sometimes face compilation errors or test failures that only manifest in specific Continuous Integration (CI) builds. 
To effectively troubleshoot these issues, it's advisable to debug directly in the target environment. The preferred method involves using Docker to locally replicate the problematic image.
However, for those without Docker access on their machines, (or for cases inherently related to the CI configuration), an alternative is to establish a connection to the CI machines. Here are the steps to do so:

Step 1: Adding a GHA to establish a connection
==============================================

First, as much as you can, try to reduce the number of jobs you're triggering by commenting out the configurations you do not require for your debugging.
Then in your branch, add the following GHA step to the `.github/build_and_test.yml` (see full documentation of the action `here <https://github.com/lhotari/action-upterm>_`).

.. code-block:: console
  - name: ssh  
      uses: lhotari/action-upterm@v1  
      with:
        ## limits ssh access and adds the ssh public key for the user which triggered the workflow
        limit-access-to-actor: true
        ## limits ssh access and adds the ssh public keys of the listed GitHub users
        limit-access-to-users: GitHubLogin

The action should be added after whichever step triggers an error. In case of a build failure it is best to add the action after the `build, test and deploy` step.
It is also important to prevent the job to exit upon failure. For instance, it is suggested to comment the following lines in the `build, test and deploy` step.

.. code-block:: console
    set -e

.. code-block:: console
    exit ${EXIT_STATUS}


You can now commit the changes and push them to your remote branch.

Step 2: Inspect the CI and grab server address
==============================================

.. code-block:: console
  Run lhotari/action-upterm@v1
  upterm
  
  Auto-generating ~/.ssh/known_hosts by attempting connection to uptermd.upterm.dev
  Pseudo-terminal will not be allocated because stdin is not a terminal.

  Warning: Permanently added 'uptermd.upterm.dev' (ED25519) to the list of known hosts.

  runner@uptermd.upterm.dev: Permission denied (publickey).

  Adding actor "GitHubLogin" to allowed users.
  Fetching SSH keys registered with GitHub profiles: GitHubLogin
  Fetched 2 ssh public keys
  Creating a new session. Connecting to upterm server ssh://uptermd.upterm.dev:22
  Created new session successfully
  Entering main loop 
  === Q16OBOFBLODJVA3TRXPL                                                                                                 
  Command:                tmux new -s upterm -x 132 -y 43                                                                 
  Force Command:          tmux attach -t upterm                                                                           

  Host:                   ssh://uptermd.upterm.dev:22                                                                     
  SSH Session:            ssh Q16oBofblOdjVa3TrXPl:ZTc4NGUxMWRiMjI5MDgudm0udXB0ZXJtLmludGVybmFsOjIyMjI=@uptermd.upterm.dev


Step 3: Connect to the machine via ssh
======================================

You can now open a terminal in your own machine and sshe to the upterm server, e.g.,


.. code-block:: console
   ssh Q16oBofblOdjVa3TrXPl:ZTc4NGUxMWRiMjI5MDgudm0udXB0ZXJtLmludGVybmFsOjIyMjI=@uptermd.upterm.dev


Step 4: Run the docker container interactively
==============================================
Once you are connected to the machine it is convenient to follow these steps to interactively run the docker container:

.. code-block:: console
   docker ps -a


The id of the existing docker container will be displayed and you can use it to commit the container.

.. code-block:: console
   docker commit <id> debug_image

and then run it interactively, e.g.

.. code-block:: console
   docker run -it --volume=/home/runner/work/GEOS/GEOS:/tmp/geos -e ENABLE_HYPRE=ON -e ENABLE_HYPRE_DEVICE=CUDA -e ENABLE_TRILINOS=OFF --cap-add=SYS_PTRACE --entrypoint /bin/bash debug_image

Step 5: Cancel the workflow
============================================== 
Once you are done, do not forget to cancel the workflow!