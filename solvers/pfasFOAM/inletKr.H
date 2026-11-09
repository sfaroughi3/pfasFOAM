{
    const scalar krInletMin = 1e-30;

    forAll(h.boundaryField(), patchi)
    {
        if (h.boundaryField().types()[patchi] == "darcyGradPressure")
        {
            krthetaf.boundaryFieldRef()[patchi] =
                max(krtheta.boundaryField()[patchi].patchInternalField(), krInletMin);
        }
    }
}
