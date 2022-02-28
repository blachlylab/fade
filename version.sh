if [ -z "${FADE_VER}" ]; then
    printf 'module _version;\nenum VERSION="%s";\n' $(git describe --tags) > source/_version.d
else
    printf 'module _version;\nenum VERSION="%s";\n' "${FADE_VER}" > source/_version.d
fi