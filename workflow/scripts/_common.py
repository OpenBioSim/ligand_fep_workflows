# Holder for common functions that need to be consistent across stages
import BioSimSpace as BSS


def runProcess(system, protocol, engine="gromacs", pmemd=False):
    """
    Given a solvated system (BSS object) and BSS protocol, run a process workflow with either
    Sander (CPU) or pmemd.cuda (GPU). NPT is typically done with GPU to save computing time.
    Returns the processed system.
    """
    engine = engine.strip().lower()

    # Create the process passing a working directory.
    if engine == "amber":
        if not pmemd:
            process = BSS.Process.Amber(system, protocol)
        elif pmemd:
            process = BSS.Process.Amber(system, protocol, is_gpu=True)

    elif engine == "gromacs":
        process = BSS.Process.Gromacs(
            system, protocol, extra_args={"-ntmpi": "1"}, ignore_warnings=True
        )

    elif engine == "openmm":
        process = BSS.Process.OpenMM(system, protocol)

    # Start the process.
    process.start()

    # Wait for the process to exit.
    process.wait()

    # Check for errors.
    if process.isError():
        print(process.stdout())
        print(process.stderr())
        raise BSS._Exceptions.ThirdPartyError("The process exited with an error!")

    # If it worked, try to get the system. No need to block, since it's already finished.
    system = process.getSystem()

    return system


if __name__ == "__main__":
    raise RuntimeError(
        "This script is not intended to be run directly. It is a module for use in other scripts."
    )
