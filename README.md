# CPNS · Tools branch

This branch collects CPNS lab tool repositories as git submodules (linked external repos). Submodules are pinned to specific commits for reproducibility, and can be updated to newer upstream commits when you choose.

Included tools (submodules):
- atcm — https://github.com/alexandershaw4/atcm
- aLogLikeFit — https://github.com/alexandershaw4/aLogLikeFit
- VariationalLaplace — https://github.com/alexandershaw4/VariationalLaplace

How to use this branch:

Clone this branch (including submodules):
```bash
git clone -b tools --recurse-submodules https://github.com/liobaberndt/CPNS.git
```

If you already cloned without submodules (or the submodule folders are empty):
```bash
git checkout tools
git submodule update --init --recursive
```

Update submodules to the latest upstream commits (submodules do not update automatically):

```bash
git checkout tools
git pull

git submodule update --remote --merge

git add .
git commit -m "Update tool submodules"
git push
```
