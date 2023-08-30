#pragma once

#ifndef PROGRESS_WINDOW_H
#define PROGRESS_WINDOW_H

#if SHOW_PROGRESS > 0

#include <CommCtrl.h>   // for progressbar

/* Returns the last Win32 error, in string format. Returns an empty string if there is no error. */
inline std::string GetLastErrorAsString()
{
    //Get the error message ID, if any.
    DWORD errorMessageID = GetLastError();
    if (errorMessageID == 0)
        return std::string(); //No error message has been recorded

    LPSTR messageBuffer = nullptr;

    //Ask Win32 to give us the string version of that message ID.
    //The parameters we pass in, tell Win32 to create the buffer that holds the message for us (because we don't yet know how long the message string will be).
    size_t size = FormatMessageA(FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
        NULL,
        errorMessageID,
        MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
        (LPSTR)&messageBuffer,
        NULL,
        NULL);

    //Copy the error message into a std::string.
    std::string message(messageBuffer, size);

    //Free the Win32's string's buffer.
    LocalFree(messageBuffer);

    return message;
};

inline LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
    workspace.progress_window = hWnd;

    switch (message)
    {
    case WM_CLOSE:
        //DestroyWindow(hWnd);
        //EndDialog(hWnd, NULL);
        break;
    case WM_DESTROY:
        //DestroyWindow(hWnd);
        //EndDialog(hWnd, NULL);
        break;
    default:
        return DefWindowProc(hWnd, message, wParam, lParam);
    }

    return 0;
};

struct ThreadArgs
{
    HINSTANCE hInstance = workspace.hInstThisDll;
    LPCWSTR lpTemplateName = MAKEINTRESOURCE(IDD_DIALOG1);
    HWND hWndParent = NULL;
    DLGPROC lpDialogFunc = (DLGPROC)WndProc;
    LPARAM dwInitParam = NULL;
};

inline DWORD WINAPI DialogBoxParamWrapper(ThreadArgs* lpParameter)
{
    if (DialogBoxParam(lpParameter->hInstance, // ÑÞÄÀ ÍÓÆÍÎ ÊËÀÑÒÜ ÒÎËÜÊÎ HINSTANCE ÝÒÎÉ DLL, ÈÍÀ×Å ÍÅ ÐÀÁÎÒÀÅÒ
        lpParameter->lpTemplateName,
        lpParameter->hWndParent,
        lpParameter->lpDialogFunc,
        lpParameter->dwInitParam) < 0)
        MessageBoxA(NULL, GetLastErrorAsString().c_str(), "Error!", MB_ICONINFORMATION | MB_OK);

    return 0;
};

#endif /* SHOW_PROGRESS > 0 */

#endif /* PROGRESS_WINDOW_H */
