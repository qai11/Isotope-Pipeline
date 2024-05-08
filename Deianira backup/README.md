# Deianira
Deianira is an attempt to fix the broken bits of MEGARA in a more generalised way, that works for all HERCULES data.

## Usage of Deianira
Make sure Fits_extension_edit.py and BArycentric_correction_iSpec.py are in the same folder. In Fits_extension_edit, replace the variable `reduced_data` with a string of the folder your .mat files from MEGARA are in. Replace `post_reduced_data` with a string of where you want to save the new .fits.

## Using git functionality
If you want to edit the functionality of the code, make a new branch.
This can be done by clicking on the branch overview on the bottom left of VSCode, or by typing
```
git checkout branch-name
```
into the terminal, inside the repo folder.

To update the branches your VSCode can access, type
```
git fetch
```
into the terminal.

To delete a branch locally that has been deleted on the server, do the following:
```
CTRL+SHIFT+P
Git: Delete branch
```
select the branch to delete, and it will be deleted locally.
To delete a branch remotely, do this on the eng-git website (here).
