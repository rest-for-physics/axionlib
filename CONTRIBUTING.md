Contributions to this library can be made directly to the `master` branch when the process, metadata, or any contributed code, is stable enough to be integrated at the next official release. This does not mean it is 100% bug free, but that the new process it is in a final state and its use would be accepted in an official data processing chain, or results will be published using that code and they need to be tagged for future reference.

Experimental developments, or developments taking longer periods of time, are encouraged to create a dedicated branch to be merged to master in the future to avoid adding to an official tag/release any unfinished job.

### Using GitHub issue tracker

New contributions are encouraged to create a new issue entry at the GitHub issue tracker in order to explain the nature of the developments, allow tracking of those developments, and invite other developers to join the discussion, evaluation and testing of the new code.

When we use the issue tracking system we must write at the commit message the issue number, so that it keeps registered at the issue tracker.

For example:

```
git commit -m "ClassName. Fixed initalization bug. Issue #1"
```

As soon as any issues remain open (and commits connected to this issue have been already added to master) those issues should be closed before fixing a new library release.

### Contribution requirements

Authors pushing new processes or metadata classes to this library will be encouraged to prompty include at least:

1. Doxygen in-code documentation describing the process pourpose and scope, including examples, and if possible, a figure ilustrating the effect of the process on event data.
2. A validation test with a minimal running test to be included at the pipeline file `.gitlab-ci.yml`. Tests will be running at https://lfna.unizar.es/iaxo/RestAxionLib.

### Fixing a new library release

In a last commit we will update manually the version found at the `CMakeLists.txt` file. I.e. from `1.0` to `1.1`,

```
set( LibraryVersion "1.1" )
```

Then we will commit and push the change,

```
git commit -m "Updating library to version 1.1"
git push
```

and we will create the new tag

```
git tag -a v1.1 -m "Fixing release 1.1"
git push --tags
```

As soon as the new release is ready, the most natural is to update the submodule at the main [framework](https://github.com/rest-for-physics/framework) in order to make official the changes in the next framework release.

If we are the main framework directory this would be achieved by doing

```
git add source/libraries/axion
git commit -m "Updating axion library submodule to version 1.1"
git push
```

### Versioning

Please, notice that the central versioning system, which guarantees code traceability, is only managed by the framework. The library version number is only used for users to identify major changes, access the release notes for the corresponding updates, and create a citable reference to be used in publications.
