#!/bin/env bash

case "$(grep "ID=" /etc/os-release)" in
        *"gentoo"*)
            echo "Gentoo system found"
        ;;
        *"debian"*)
            echo "Debian system found"
        ;;
        *"arch"*)
            echo "Arch system found"
        ;;
        *"ubuntu"*)
            echo "Ubuntu system found"
        ;;
        *)
            echo $(grep "ID=" /etc/os-release) "system found: trying defaults"
        ;;
esac
