# Development Process

To date, we've had only one developer and the process has simply been:
- Push tested commits to the development branch
- At appropriate points, create new releases using the steps described below
- At each release point move all commits to master via fast-forward rebase
  and push.  Thus the latest commit on master is always a release point.

Plan is to move to a pull request model if we get more developers.

## Release Numbering

Releases are numbered using three components: X.Y.Z.
- Increment Z for bug fixes and minor changes that don't change the API
  or intended semantics in any way.
- Increment Y and set Z to zero when new functions or new behavior is
  added.
- For very significant new but still backward compatible functionality
  (rare), increment X and set Y and Z to zero.

Note that we MUST NOT ever remove functionality or otherwise break binary
or source-level compatibility, even when X is incremented.  The reason
is that X.Y.Z is the package version, thus upgrade to any higher version
should safely be able to remove the previous version without breaking
existing programs.

If a change cannot avoid breaking compatiblity, it must be released as
a whole new library, e.g. libvna2-1.0.0.  This approach provides maximum
coexistence between incompatible library versions.

Exception: we have broken compatiblity a few times in the early 0.y.z
releases.  At this point, however, we will not break compatibliity.

TODO: add in library versioning rules

## Creating a new release

- Run make -j12 check
- Run make distcheck
- Run make rpm
- Push to origin/development
- Make sure CI/CD build succeeds
- On the development branch, update the version in the AC_INIT macro
  of configure.ac.  
- Commit with a message in this format "version 0.0.3"
- Create an annotated tag.  Example:

```
git tag -a v0.0.3 -m "v0.0.3"
```

- Checkout the master branch, add the release and push to github

```
git checkout master
git rebase development
git push
git push --tags
```
